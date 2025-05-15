
"""
    on_grid(plasma::T, p0::T) where {T<:Real}

Checks if p0 is on the plasma grid
"""
function on_grid(plasma::Plasma, p0::Vector{T}) where {T<:Real}
    R0 = hypot(p0[1], p0[2])
    z0 = p0[3]
    return !(plasma.R_coords[1] <= R0 <= plasma.R_coords[end] && plasma.Z_coords[1] <= z0 <= plasma.Z_coords[end])
end

"""
    first_point(p0::T, k0::T, plasma::T)

Finds the first point on the plasma grid along the ray
"""
function first_point(plasma::Plasma, p0::Vector{T}, N0::Vector{T}) where {T<:Real}
    if on_grid(plasma, p0)
        p_plasma = IMAS.ray_rect_torus_intersect(p0, N0, 
                                                 [plasma.R_coords[1], plasma.R_coords[end]],
                                                 [plasma.Z_coords[1], plasma.Z_coords[end]])[1]
    else
        p_plasma = p0
    end
    return p_plasma
end

function refraction_equations!(F::Vector{T}, x::Vector{T}, X::T, Y::T, n::Vector{T}, b::Vector{T}, mode::Integer)  where {T<:Real}
    Nz = sqrt(1 - hypot(x[2:4])^2)
    N = [x[2], x[3], Nz]
    N_par = dot_product(N, b) * x[1]
    n_dot_N = -dot_product(n, N)
    F[1] = refractive_index_sq(X, Y, N_par, mode)
    F[2:5] = N/x[1] + (1/x[1] * n_dot_N -  sqrt(1.0 - 1/x[1]^2 * (1.0 - n_dot_N^2))) * n
    F[1] = F[1] - x[1]
    F[2:5] = F[2:5] - N
end

function vacuum_plasma_refraction(plasma::Plasma, p_plasma::T, N0::T, omega::Real, mode::Integer) where {T<:Real}
    X, Y, N_par = eval_plasma(plasma, p_plasma[1], p_plasma[2], p_plasma[3], N0[1], N0[2], N0[3], omega)
    # Estimate refractive index for perpendicular propagation
    N_est = refractive_index_sq(X, Y, 0.0, mode)
    # Use this estimate to check if we get an immediate reflection at vaccum -> plasma boundary
    if imag(N_est) != 0 || N_est < 0
        return false, nothing
    end
    # Calculate flux surface normal
    R = hypot(p_plasma[1], p_plasma[2])
    phi = atan(p_plasma[2], p_plasma[1])
    dpsi_dR, dpsi_dz = Interpolations.gradient(plasma.psi_norm_spline, R, p_plasma[3])
    n = Vector{Float64}(3)
    b = B_spline(plasma, R, phi, p_plasma[3])
    b = b / hypot(b)
    # This assumes dpsi_dphi is zero!
    n[1] = dpsi_dR * p_plasma[1] / R
    n[2] = dpsi_dR * p_plasma[2] / R
    n[3] = dpsi_dz
    n = n / hypot(n)
    x0 = [N_est, N0[1], N0[2]]
    x = nlsolve((F,x) ->f!(F, x, X, Y, n, b, mode), x0, autodiff = :forward)
    n_plasma = x[1] * [x[2], x[3], sqrt(1 - hypot(x[2:4])^2)]
    return n_plasma
end
"""
    solve(plasma::Plasma, r0::T, ϕ0::T, z0::T, nϕ0::T, θ_injection::T, freq::T, tmax::Float64) where {T<:Real}

Launch the ray in a given plasma, given some initial conditions
"""
function make_ray(plasma::Plasma, x0::T, y0::T, z0::T, steering_angle_tor::T, steering_angle_pol::T, freq::T, mode::Integer, s_max::Float64) where {T<:Real}
    N_vacuum = collect(IMAS.pol_tor_angles_2_vector(steering_angle_tor, steering_angle_pol))
    
    p_plasma = first_point(plasma, [x0, y0, z0], N_vacuum)
    n_plasma = vacuum_plasma_refraction(plasma, p_plasma, N_vacuum, ω, mode)
    
    @assert abs(λ(p_plasma, n_plasma)) < 1e-12 "Initial state is not on λ = 0 surface"

    u0 = [p_plasma; n_plasma] # concatenate for the solver
    sspan = (0.0, 10.0)  # arc‑length interval

    ######################################################################
    # 1.  Symbols                                                         #
    ######################################################################

    @parameters s, ω, p, m               # independent variable “arc length”
    @variables  (x(s))[1:3] (N(s))[1:3]  # x(s) = (x₁,x₂,x₃)ᵀ , N(s) = (N₁,N₂,N₃)ᵀ
    params = [ω => 2 * pi *freq, p => plasma, m => mode]
    D = Differential(s)

    ######################################################################
    # 2.  The level‑set  λ(x,N)                                           #
    ######################################################################
    # <<<  EDIT THIS FUNCTION TO MATCH YOUR PROBLEM  >>>
    λ(x,N) = (x,N) -> dispersion_relation(x, N, p, ω, mode) 

    ######################################################################
    # 3.  Symbolic gradients ∂λ/∂x  and  ∂λ/∂N                            #
    ######################################################################
    λ_expr    = λ(x, N)                             # symbolic expression
    ∇xλ       = expand_derivatives.(Symbolics.gradient(λ_expr, x))   # length‑3 Vector{Num}
    ∇Nλ       = expand_derivatives.(Symbolics.gradient(λ_expr, N))

    norm_∇Nλ  = sqrt(sum(∇Nλ .^ 2))                # ‖∂λ/∂N‖

    ######################################################################
    # 4.  Differential equations                                          #
    ######################################################################
    eqs = [
        [D(x[i]) ~  ∇Nλ[i] / norm_∇Nλ                for i in 1:3],
        [D(N[i]) ~ -∇xλ[i] / norm_∇Nλ                for i in 1:3]
    ]

    hamiltonian = ODESystem(eqs, t = s, name = :Characteristics)

    prob = ODEProblem(structural_simplify(hamiltonian), u0, sspan, params)
    sol  = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-10)
end

# function out_of_bounds(u, t, integrator)
#     r = u[1]
#     z = u[2]
#     return r < plasma.R_coords[1] || r > plasma.R_coords[end] || z < plasma.Z_coords[1] || z > plasma.Z_coords[end]
# end