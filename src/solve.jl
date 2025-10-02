


"""
    Base.@kwdef struct TorJParams

Provides a structure for code parameters. All values have defaults.
    - s_max: # Stopping condition for the raytracer if there is insufficient absorption [m]. Defaults to 1 m
    - N_rings: Number of radial rings of rays in the beam. Defaults to 3
    - min_azimuthal_points: Number of poloidal rays in the innermost rings. The outer rings will have more than this. Defaults to 5
"""

Base.@kwdef struct TorJParams
    s_max:: Float64 = 1.0 
    N_rings:: Int64 = 3 
    min_azimuthal_points:: Int64 = 5
    N_dP_dV::Int64 = 1000
end



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
    p_plasma += Roots.find_zero(g, (0.0, 0.5), Roots.Bisection(); xtol=1e-6).* N0
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
    solver_results = NLsolve.nlsolve((F,N) -> refraction_equations!(F, N, X, Y, N0 / LinearAlgebra.norm(N0), n, b, mode), 
                             x0, autodiff = :forward; ftol = 5.e-13)
    return solver_results.zero
end
"""
Simple routine that enforces non-negativity for the power
"""
function affect!(integrator)
    # Project only specific components to zero
    if integrator.u[7] < 0
        integrator.u[7] = 0.0
    end
end

function gradΛ!(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, plasma::Plasma, ω::Real, mode::Integer)
    # Struggling to get rid of this allocation.
    # It should be easy to use `du` as a temporary array but 
    # in practice it always leads to errors
    ForwardDiff.gradient!(@view(du[4:6]), x -> dispersion_relation(x, u[4:6], plasma, ω, mode), u[1:3])
    ForwardDiff.gradient!(@view(du[1:3]), N -> dispersion_relation(u[1:3], N, plasma, ω, mode), u[4:6])
    norm_∂Λ_∂N = LinearAlgebra.norm(du[1:3])
    du[1:6] ./= norm_∂Λ_∂N
    du[4:6] .= -du[4:6]
    du[7] = -u[7] * α_approx(u[1:3], u[4:6], plasma, ω, mode)
end    

"""
    sys!(du, u, p, s, plasma::Plasma, ω::Real, mode::Integer)

System of ODEs for ray tracing with absorption.
Defines du/ds = f(u, s) where u = [x, N, P] contains position, refractive index vector, and power.

# Arguments
- `du`: Derivative vector [dx/ds, dN/ds, dP/ds]
- `u`: State vector [x, N, P] (position, refractive index, power)
- `p`: Parameters (unused)
- `s`: Arclength parameter
- `plasma`: Plasma configuration
- `ω`: Angular frequency
- `mode`: Polarization mode (+1 for X-mode, -1 for O-mode)
"""
function sys!(du, u, p, s, plasma::Plasma, ω::Real, mode::Integer)
    gradΛ!(du, u, plasma, ω, mode)
end

"""
    make_ray(plasma::Plasma, x0::AbstractVector{<:T}, N_vacuum::AbstractVector{<:T}, f::T, mode::Integer, s_max::Float64) where {T<:Real}

Launch a single ray through the plasma and compute its trajectory with absorption.

# Arguments
- `plasma`: Plasma configuration containing density, temperature, and magnetic field
- `x0`: Initial position in vacuum [x, y, z] (meters)
- `N_vacuum`: Initial refractive index vector in vacuum
- `f`: Wave frequency (Hz)
- `mode`: Polarization mode (+1 for X-mode, -1 for O-mode)
- `s_max`: Maximum integration distance (meters)
- `psi_dP_dV`: Psi axis for power deposition profile

# Returns
- `s`: Arclength array along the trajectory
- `u`: Position trajectory as vector of 3D points
- `P_beam`: Beam power along the trajectory
"""
function make_ray(plasma::Plasma, x0::AbstractVector{<:T}, N_vacuum::AbstractVector{<:T}, 
                  f::T, mode::Integer, s_max::Float64, psi_dP_dV::AbstractVector{<:T}) where {T<:Real}
    p_plasma = first_point(plasma, x0, N_vacuum)
    @assert evaluate(plasma.psi_norm_spline, p_plasma) <= plasma.psi_prof_max

    N_plasma = vacuum_plasma_refraction(plasma, p_plasma, N_vacuum, 2 * pi *f, mode)
    initial_lambda = abs(dispersion_relation(p_plasma, N_plasma, plasma, 2 * pi *f, mode))
    @assert initial_lambda < 1e-12 "Initial state is not on λ = 0 surface. λ=$initial_lambda"

    # Initial condition: [position, refractive_index, power]
    u0 = [p_plasma; N_plasma; 1.0]
    N_steps = 100
    s_step = s_max / Float64(N_steps)
    # We will add the arc_length for the step in vacuum in the end
    s0 = LinearAlgebra.norm(p_plasma .- x0)
    s = vec([0.0, s0])
    P_beam = vec([1.0, 1.0])
    dP_ds = vec([0.0, 0.0])
    # Add the first two points
    u = vec([vec(x0), p_plasma])
    for i in 1:N_steps
        prob = DifferentialEquations.ODEProblem((du, u, p, s) -> sys!(du, u, p, s, plasma, 2.0 * pi *f, mode), u0, 
                                                (Float64(i-1)*s_step + s0, Float64(i)*s_step + s0), 
                                                DifferentialEquations.OwrenZen3(); 
                                                dtmax=1.e-4, abstol=1.e-6, reltol=1.e-6)
        # Ensure that the power does not oscillate around 0.0
        condition(u, t, integrator) = u[7] < 0
        cb = DifferentialEquations.ContinuousCallback(condition, affect!)
        sol = DifferentialEquations.solve(prob, callback=cb)
        u0 = sol.u[end]
        append!(s, sol.t[2:end])
        # Extract position vectors (first 3 components) from each solution step
        # Discard the first step which is a repetition
        for solution_step in sol.u[2:end]
            push!(u, solution_step[1:3])  # Only store position, not direction
            push!(P_beam, solution_step[7])
            # Store the dP/ds needed for power deposition.
            # It is difficult to reuse information from the ODE integration because `sol.u` is an interpolation
            push!(dP_ds, solution_step[7] * α_approx(solution_step[1:3], solution_step[4:6], plasma, 2.0 * pi *f, mode))
        end
        # Check if we are still in the plasma otherwise stop propagation
        if evaluate(plasma.psi_norm_spline, u0[1:3]) > 1.0 break end
        # Check if the beam has any power left, otherwise stop propagation
        if sol.u[end][7] < 1.e-6 break end
    end
    # Finally add the offset from the vacuum step
    dP_dV_ray, deposited_power = power_deposition_profile(plasma, s, u, dP_ds, psi_dP_dV)
    return s, u, P_beam, dP_dV_ray, deposited_power
end

"""
    spawn_beam_tasks(plasma, r, phi, z, steering_angle_tor, steering_angle_pol, spot_size, inverse_curvature_radius, f, mode, s_max, psi_dP_dV; kwargs...)

Spawn Dagger tasks for multiple rays forming a beam.

# Arguments
- `plasma`: Plasma configuration containing density, temperature, and magnetic field
- `r`: Radial distance of beam center (meters)
- `phi`: Toroidal angle of beam center (radians)
- `z`: Vertical position of beam center (meters)
- `steering_angle_tor`: Toroidal steering angle (radians)
- `steering_angle_pol`: Poloidal steering angle (radians)
- `spot_size`: Beam spot size parameter (meters)
- `inverse_curvature_radius`: Inverse radius of curvature (1/meters)
- `f`: Wave frequency (Hz)
- `mode`: Polarization mode (+1 for X-mode, -1 for O-mode)
- `s_max`: Maximum integration distance (meters)
- `psi_dP_dV`: Psi axis for power deposition profile

# Returns
- `ray_tasks`: Vector of Dagger task objects
- `ray_weights`: Ray weights for beam integration
"""
function spawn_beam_tasks(plasma::Plasma, r::T, phi::T, z::T, steering_angle_tor::T, steering_angle_pol::T, spot_size::T, 
                          inverse_curvature_radius::T, f::T, mode::Integer, s_max::Float64, psi_dP_dV::Vector{T}; kwargs...) where {T<:Real}
    N0 = collect(IMAS.pol_tor_angles_2_vector(steering_angle_pol, steering_angle_tor))
    x0 = zeros(Float64,3)
    x0[1] = r * cos(phi)
    x0[2] = r * sin(phi)
    x0[3] = z
    ray_positions, ray_directions, ray_weights = TorJ.launch_peripheral_rays(
        x0, N0, spot_size, inverse_curvature_radius, f; kwargs...)
    ray_tasks = []
    Dagger.spawn_bulk() do
        for i in eachindex(ray_weights)
            push!(ray_tasks, Dagger.@spawn make_ray(plasma, ray_positions[i,:], ray_directions[i,:], 
                                                    f, mode, s_max, psi_dP_dV) )
        end
    end
    return ray_tasks, ray_weights
end

"""
    fetch_beam_results(ray_tasks, ray_weights, psi_dP_dV)

Fetch and process results from spawned beam tasks.

# Arguments
- `ray_tasks`: Vector of Dagger task objects
- `ray_weights`: Ray weights for beam integration
- `psi_dP_dV`: Psi axis for power deposition profile

# Returns
- `arc_lengths`: Vector of arclength arrays (one per ray)
- `trajectories`: Vector of trajectory arrays (one per ray)  
- `ray_powers`: Vector of beam power arrays (one per ray)
- `dP_dV`: Integrated power deposition profile
- `deposited_power`: Total deposited power
- `ray_weights`: Ray weights for beam integration
"""
function fetch_beam_results(ray_tasks, ray_weights, psi_dP_dV)
    ray_results = []
     # Collect results preserving order (ensures correspondence with ray_weights)
    for i in eachindex(ray_tasks)
        push!(ray_results, Dagger.fetch(ray_tasks[i]))
    end
    # Separate the three return values
    arc_lengths = [result[1] for result in ray_results]
    trajectories = [result[2] for result in ray_results]
    ray_powers = [result[3] for result in ray_results]
    dP_dV_ray  = [result[4] for result in ray_results]
    deposited_power_ray  = [result[5] for result in ray_results]
    
    dP_dV = zeros(Float64, length(psi_dP_dV))
    # From the power deposition profile
    deposited_power = 0.0
    # From the rays directly
    for i_ray in eachindex(ray_weights)
        dP_dV .+= dP_dV_ray[i_ray] * ray_weights[i_ray]
        deposited_power += deposited_power_ray[i_ray] * ray_weights[i_ray]
    end
    return arc_lengths, trajectories, ray_powers, dP_dV, deposited_power, ray_weights 
end

"""
    make_beam(plasma, r, phi, z, steering_angle_tor, steering_angle_pol, spot_size, inverse_curvature_radius, f, mode, s_max)

Launch multiple rays forming a beam and compute their trajectories with absorption in parallel using Dagger.

# Arguments
- `plasma`: Plasma configuration containing density, temperature, and magnetic field
- `r`: Radial distance of beam center (meters)
- `phi`: Toroidal angle of beam center (radians)
- `z`: Vertical position of beam center (meters)
- `steering_angle_tor`: Toroidal steering angle (radians)
- `steering_angle_pol`: Poloidal steering angle (radians)
- `spot_size`: Beam spot size parameter (meters)
- `inverse_curvature_radius`: Inverse radius of curvature (1/meters)
- `f`: Wave frequency (Hz)
- `mode`: Polarization mode (+1 for X-mode, -1 for O-mode)
- `s_max`: Maximum integration distance (meters)
- `psi_dP_dV`: Psi axis for power deposition profile

# Returns
- `arc_lengths`: Vector of arclength arrays (one per ray)
- `trajectories`: Vector of trajectory arrays (one per ray)  
- `ray_powers`: Vector of beam power arrays (one per ray)
- `dP_ds`: Vector of change of power as function of s (one per ray)
- `ray_weights`: Ray weights for beam integration
"""
function make_beam(plasma::Plasma, r::T, phi::T, z::T, steering_angle_tor::T, steering_angle_pol::T, spot_size::T, 
                   inverse_curvature_radius::T, f::T, mode::Integer, s_max::Float64, psi_dP_dV::Vector{T}; kwargs...) where {T<:Real}
    ray_tasks, ray_weights = spawn_beam_tasks(plasma, r, phi, z, steering_angle_tor, steering_angle_pol, spot_size, 
                                              inverse_curvature_radius, f, mode, s_max, psi_dP_dV; kwargs...)
    return fetch_beam_results(ray_tasks, ray_weights, psi_dP_dV)
end

function process_ids!(dd::IMAS.dd; torj_params::TorJParams = TorJParams())
    nbeam = length(dd.ec_launchers.beam)
    if nbeam < 1
        return nbeam
    end

    plasma = Plasma(dd)
    
    psi_norm_dP_dV = Vector(LinRange(0.0, 1.0, torj_params.N_dP_dV))
    beam_tasks = []
    beam_weights = []
    # This spawns the beam tasks serially so that we can parallelize over rays first
    Dagger.spawn_sequential() do
        for i_beam in 1:nbeam
            beam = dd.ec_launchers.beam[i_beam]
            f = IMAS.@ddtime(beam.frequency.data)

            # Elliptical beam parameters (not supported yet) -> circular beam
            inverse_curvature_radius = 0.5*(beam.phase.curvature[1, 1] + beam.phase.curvature[2, 1])
            spot_size = 0.5*(beam.spot.size[1, 1] + beam.spot.size[2, 1])
            
            # Sign convention opposite to IMAS
            mode = -beam.mode
            ray_tasks, ray_weights = spawn_beam_tasks(plasma, 
                    beam.launching_position.r[1], beam.launching_position.phi[1], 
                    beam.launching_position.z[1], 
                    IMAS.@ddtime(dd.ec_launchers.beam[i_beam].steering_angle_tor), 
                    IMAS.@ddtime(dd.ec_launchers.beam[i_beam].steering_angle_pol), 
                    spot_size, inverse_curvature_radius, f, mode, 
                    torj_params.s_max, psi_norm_dP_dV; N_rings = torj_params.N_rings, 
                    min_azimuthal_points = torj_params.min_azimuthal_points)
            push!(beam_tasks, ray_tasks)
            push!(beam_weights, ray_weights)
        end
    end

    # Fetch results from all beams
    beam_results = []
    for i_beam in eachindex(beam_tasks)
        result = fetch_beam_results(beam_tasks[i_beam], beam_weights[i_beam], psi_norm_dP_dV)
        push!(beam_results, result)
    end

    # Prepare the output
    
    resize!(dd.waves.coherent_wave, nbeam)
    for i_beam in eachindex(beam_results)
        arc_lengths, trajectories, ray_powers, dP_dV, deposited_power, ray_weights = beam_results[i_beam]
        beam = dd.ec_launchers.beam[i_beam]
        ps_beam = dd.pulse_schedule.ec.beam[i_beam]
        power_launched = IMAS.@ddtime(ps_beam.power_launched.reference)

        wv = dd.waves.coherent_wave[i_beam]
        wv.identifier.antenna_name = beam.name
        wv.identifier.type.description = "TorJ"
        wv.identifier.type.name = "EC"
        wv.identifier.type.index = 1
        wv.wave_solver_type.index = 1 # BEAM/RAY TRACING

        wvg = resize!(wv.global_quantities) # global_time
        wvg.frequency = IMAS.@ddtime(beam.frequency.data)
        # torj never sees the launched power so we handle it here
        wvg.electrons.power_thermal = deposited_power * power_launched
        wvg.power = deposited_power * power_launched

        wv1d = resize!(wv.profiles_1d) # global_time
        wv.profiles_1d[1].time = IMAS.@ddtime(dd.equilibrium.time)
        wv1d.grid.rho_tor_norm = plasma.rho_tor_psi_norm_spline(psi_norm_dP_dV)
        wv1d.grid.psi = psi_norm_dP_dV
        wv1d.power_density = dP_dV .* power_launched

        source = resize!(dd.core_sources.source, :ec, "identifier.name" => beam.name; wipe=false)
        IMAS.new_source(
            source,
            source.identifier.index,
            beam.name,
            wv1d.grid.rho_tor_norm,
            wv1d.grid.volume,
            wv1d.grid.area;
            electrons_energy=wv1d.power_density)

        # LOOP OVER RAYS
        wvb = resize!(wv.beam_tracing) # global_time
        resize!(wvb.beam, length(ray_weights)) # Five beams/per gyrotron
        if power_launched > 0.0
            for i_ray in eachindex(ray_weights)
                wvb.beam[i_ray].power_initial = ray_weights[i_ray] * power_launched
                s = arc_lengths[i_ray]
                x = [point[1] for point in trajectories[i_ray]]
                y = [point[2] for point in trajectories[i_ray]]
                z = [point[3] for point in trajectories[i_ray]]
                r = hypot.(x, y)
                wvb.beam[i_ray].length = s
                wvb.beam[i_ray].position.r = r
                wvb.beam[i_ray].position.phi = atan.(y,x)
                wvb.beam[i_ray].position.z = z
                wvb.beam[i_ray].electrons.power = (1.0 .- ray_powers[i_ray]) .* (ray_weights[i_ray] * power_launched)
            end
        end

    end # LOOP OVER BEAMS (LAUNCHERS)

    IMAS.@ddtime(dd.waves.code.output_flag = 0) # NO ERROR
    return dd
end