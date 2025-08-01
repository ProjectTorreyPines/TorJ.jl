
"""
    on_grid(plasma::T, p0::T) where {T<:Real}

Checks if p0 is on the plasma grid
"""
function on_grid(plasma::Plasma, p0::AbstractVector{T}) where {T<:Real}
    R0 = hypot(p0[1], p0[2])
    z0 = p0[3]
    return plasma.R_coords[1] <= R0 <= plasma.R_coords[end] && plasma.Z_coords[1] <= z0 <= plasma.Z_coords[end]
end

"""
    first_point(p0::T, k0::T, plasma::T)

Finds the first point on the plasma grid with plasma profile data along the ray
"""
function first_point(plasma::Plasma, p0::AbstractVector{T}, N0::Vector{T}) where {T<:Real}
    if on_grid(plasma, p0)
        p_plasma = p0
    else
        t = IMAS.toroidal_intersection([plasma.R_coords[1], plasma.R_coords[end], plasma.R_coords[end], plasma.R_coords[1], plasma.R_coords[1]], 
                                       [plasma.Z_coords[1], plasma.Z_coords[1], plasma.Z_coords[end], plasma.Z_coords[end], plasma.Z_coords[1]], 
                                        p0, N0)
        p_plasma = p0 .+ N0 .* t
    end
    g(t) = evaluate(plasma.psi_norm_spline, p_plasma .+ t .* N0) - plasma.psi_prof_max
    # TODO find better estimate for upper limit of distance between first point on grid and plasma.psi_prof_max
    p_plasma += find_zero(g, (0.0, 0.5), Bisection(); xtol=1e-6).* N0
    psi_ref = evaluate(plasma.psi_norm_spline, p_plasma)
    # If root finding failed horribly we catch it here
    @assert abs(psi_ref - plasma.psi_prof_max) < 1.e-6
    if psi_ref > plasma.psi_prof_max
        # Move a tiny bit in
        p_plasma += 2.0 * (psi_ref - plasma.psi_prof_max) .* N0
    end
    return p_plasma
end

function refraction_equations!(F::AbstractVector{<:Real}, N::AbstractVector{<:Real}, X::T, Y::T, n0::AbstractVector{T}, 
                               n::AbstractVector{T}, b::AbstractVector{T}, mode::Integer) where {T<:Real}
    n_dot_N = -LinearAlgebra.dot(n, n0)
    # Non-linearity here because we use N
    N_par = sum(N .* b)
    Ns = sqrt(refractive_index_sq(X, Y, N_par, mode))
    F .= n0 ./ (Ns) .+ (1.0/Ns * n_dot_N .-  sqrt(1.0 - 1.0/Ns^2 * (1.0 - n_dot_N^2))) .* n # This should return a unit vector
    F .= F .* F .* Ns^2 # Multiply with Ns and square so that the absolute value represents the plasma refractive index sqaured
    F .-= N .* N
end

function vacuum_plasma_refraction(plasma::Plasma, p_plasma::AbstractVector{T}, 
                                  N0::AbstractVector{T}, omega::T, mode::Integer) where {T<:Real}
    X, Y, N_par, b = eval_plasma(plasma, p_plasma, N0, omega)
    # Estimate refractive index for perpendicular propagation
    N_est = refractive_index_sq(X, Y, 0.0, mode)
    # Use this estimate to check if we get an immediate reflection at vaccum -> plasma boundary
    if N_est <= 0
        return false, nothing
    end
    N_est = sqrt(N_est)
    # Calculate flux surface normal
    R = hypot(p_plasma[1], p_plasma[2])
    dpsi_dR, dpsi_dz = gradient(plasma.psi_norm_spline, R, p_plasma[3])
    n = zeros(Float64, 3)
    # This assumes dpsi_dphi is zero!
    n[1] = dpsi_dR * p_plasma[1] / R
    n[2] = dpsi_dR * p_plasma[2] / R
    n[3] = dpsi_dz
    n = n ./ LinearAlgebra.norm(n)
    x0 = N0 .* N_est
    solver_results = nlsolve((F,N) -> refraction_equations!(F, N, X, Y, N0 / LinearAlgebra.norm(N0), n, b, mode), 
                             x0, autodiff = :forward; ftol = 1.e-12)
    return solver_results.zero
end

function gradΛ!(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, plasma::Plasma, ω::Real, mode::Integer)
    # Struggling to get rid of this allocation.
    # It should be easy to use `du` as a temporary array but 
    # in practice it always leads to errors
    ForwardDiff.gradient!(@view(du[4:6]), x -> dispersion_relation(x, u[4:6], plasma, ω, mode), u[1:3])
    ForwardDiff.gradient!(@view(du[1:3]), N -> dispersion_relation(u[1:3], N, plasma, ω, mode), u[4:6])
    norm_∂Λ_∂N = LinearAlgebra.norm(du[1:3])
    du ./= norm_∂Λ_∂N
    du[4:6] .= -du[4:6]
end    

# Define the system of ODEs: du/ds = f(u, s)
function sys!(du, u, p, s, plasma::Plasma, ω::Real, mode::Integer)
    gradΛ!(du, u, plasma, ω, mode)
end

"""
    make_ray(plasma::Plasma, x0::AbstractVector{<:T}, N_vacuum::::AbstractVector{<:T}, freq::T, mode::Integer, s_max::Float64) where {T<:Real}

Launch the ray in a given Plasma and ray trajectory.
mode: +1 X-mode, -1 O-mode
"""
function make_ray(plasma::Plasma, x0::AbstractVector{<:T}, N_vacuum::AbstractVector{<:T}, freq::T, mode::Integer, s_max::Float64) where {T<:Real}
    # println(N_vacuum)
    p_plasma = first_point(plasma, x0, N_vacuum)
    # Make sure we are on the grid and within the bounds of the profiles for the first step
    @assert evaluate(plasma.psi_norm_spline, p_plasma) <= plasma.psi_prof_max

    N_plasma = vacuum_plasma_refraction(plasma, p_plasma, N_vacuum, 2 * pi *freq, mode)
    # println(N_plasma)
    @assert abs(dispersion_relation(p_plasma, N_plasma, plasma, 2 * pi *freq, mode)) < 1e-12 "Initial state is not on λ = 0 surface"

    u0 = [p_plasma; N_plasma] # concatenate for the solver
    # Compute ∂Λ/∂x and ∂Λ/∂N

    # Initial condition u₀ = [x₀, N₀] that satisfies Λ(x₀, N₀) ≈ 0
    # Solve the ODE
    N_steps = 100
    s_step = s_max / Float64(N_steps)
    # Add the first two points
    u = vec([vec(x0), p_plasma])
    for i in 1:N_steps
        prob = ODEProblem((du, u, p, s) -> sys!(du, u, p, s, plasma, 2.0 * pi *freq, mode), u0, 
                          (Float64(i-1)*s_step, Float64(i)*s_step), OwrenZen3(); 
                          dtmax=1.e-4, abstol=1.e-6, reltol=1.e-6)
        sol = solve(prob)
        u0 = sol.u[end]
        append!(u, sol.u)
        if evaluate(plasma.psi_norm_spline, u0[1:3]) > 1.0 break end
        # X, Y, N_par, b = eval_plasma(plasma, u0[1:3], u0[4:6], 2.0 * pi *freq)
        # R = hypot(u0[1], u0[2])
        # dX_dx_ana_test = dX_dx(plasma, u0[1:3], 2.0 * pi *freq)
        # dX_dx_test = ForwardDiff.gradient((x) -> wrap_eval_plasma(plasma, x, u0[4:6], 2.0 * pi *freq, 1), u0[1:3])
        # println("Analytical dX_dx: ", dX_dx_ana_test)
        # println("Forward diff dX_dx: ", dX_dx_test)
        # println("Analytical - forward diff dX_dx: ", dX_dx_ana_test - dX_dx_test)
        # println("Normalized error: ", 0.5*(dX_dx_ana_test .- dX_dx_test)./(dX_dx_ana_test .+ dX_dx_test))
        # N_s = sqrt(refractive_index_sq(X, Y, N_par, mode))
        # error = dispersion_relation(u0, plasma, 2.0 * pi *freq, mode)
        # println("$R, $X, $Y, $N_par, $N_s, $error")
    end
    return u
end


function make_beam(plasma::Plasma, r::T, phi::T, z::T, steering_angle_tor::T, steering_angle_pol::T, spot_size::T, 
                   inverse_curvature_radius::T, freq::T, mode::Integer, s_max::Float64) where {T<:Real}
    N0 = collect(IMAS.pol_tor_angles_2_vector(steering_angle_pol, steering_angle_tor))
    x0 = zeros(Float64,3)
    x0[1] = r * cos(phi)
    x0[2] = r * sin(phi)
    x0[3] = z
    ray_positions, ray_directions, ray_weights = TorJ.launch_peripheral_rays(
        x0, N0, spot_size, 3, inverse_curvature_radius, freq)
    trajectories = Vector{Vector{Vector{Float64}}}()
    for i in 1:length(ray_weights)
        u = make_ray(plasma, ray_positions[i,:], ray_directions[i,:], freq, mode, s_max)
        push!(trajectories, u)
    end
    return trajectories, ray_weights 
end
