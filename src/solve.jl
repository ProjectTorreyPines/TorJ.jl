using IMAS

function get_first_point(x0::T, y0::T, z0::T, steering_angle_tor::T, steering_angle_pol::T)
    kx, ky, kz = IMAS.pol_tor_angles_2_vector(steering_angle_tor, steering_angle_pol)
    
end
"""
    solve(plasma::Plasma, r0::T, ϕ0::T, z0::T, nϕ0::T, θ_injection::T, freq::T, tmax::Float64) where {T<:Real}

Launch the ray in a given plasma, given some initial conditions
"""
function solve(plasma::Plasma, x0::T, y0::T, z0::T, steering_angle_tor::T, steering_angle_pol::T, freq::T, s_max::Float64) where {T<:Real}
    ω0 = 2π * freq

    # Define the parameters dictionary, if you have any parameters
    params = [ω => ω0, p => plasma]
    kx = -cos(pol_angle) * cos(tor_angle)
    ky = -sin(tor_angle)
    kz = sin(pol_angle) * cos(tor_angle)
    # figure out kr0 and kz0 based on injection angle
    res = Optim.optimize(kz -> dispersion_relation(plasma, x0, y0, z0, 0.0, 0.0, kz, ω0), 0.0, 1E6, Optim.GoldenSection(), rel_tol=1E-3)
    k0 = res.minimizer
    kϕ0 = nϕ0 / r0
    kpol0 = sqrt(k0^2 - kϕ0^2)
    kr0 = kpol0 * cos(θ_injection)
    kz0 = kpol0 * sin(θ_injection)

    # Define initial conditions and parameter values
    ics = [
        r => r0,
        ϕ => ϕ0,
        z => z0,
        kr => kr0,
        nϕ => nϕ0,
        kz => kz0,
        t => 0.0
    ]

    # Set up the problem
    tspan = (0.0, tmax)
    prob = ODEProblem(ray_model, ics, tspan, params)

    # Solve the problem using a suitable solver from OrdinaryDiffEq
    return ModelingToolkit.solve(prob, Tsit5(), abstol=1e-6, reltol=1e-6)
end

# function out_of_bounds(u, t, integrator)
#     r = u[1]
#     z = u[2]
#     return r < plasma.R_coords[1] || r > plasma.R_coords[end] || z < plasma.Z_coords[1] || z > plasma.Z_coords[end]
# end