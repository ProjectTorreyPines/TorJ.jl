module TorJ

import IMAS
import Optim
import Interpolations: CubicSplineInterpolation, Flat, AbstractInterpolation
import ForwardDiff: derivative
using ModelingToolkit, OrdinaryDiffEq

struct Plasma{R<:AbstractRange,I<:AbstractInterpolation}
    R_coords::R
    Z_coords::R
    ne_spline::I
    Br_spline::I
    Bz_spline::I
    Bϕ_spline::I
end

function B_spline(plasma::Plasma, r, z)
    Br = plasma.Br_spline(r, z)
    Bϕ = plasma.Bϕ_spline(r, z)
    Bz = plasma.Bz_spline(r, z)
    B = sqrt(Br^2 + Bϕ^2 + Bz^2)
    return B
end

function Plasma(dd::IMAS.dd{T}) where {T}
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    eqt2d = eqt.profiles_2d[1]
    ne_data = 10.0 .^ IMAS.interp1d(cp1d.grid.psi, log10.(cp1d.electrons.density)).(eqt2d.psi)
    Br_data = eqt2d.b_field_r
    Bz_data = eqt2d.b_field_z
    Bϕ_data = eqt2d.b_field_tor

    R_coords = LinRange(eqt2d.grid.dim1[1], eqt2d.grid.dim1[end], length(eqt2d.grid.dim1))
    Z_coords = LinRange(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end], length(eqt2d.grid.dim2))

    # Interpolation objects
    ne_spline = CubicSplineInterpolation((R_coords, Z_coords), ne_data; extrapolation_bc=Flat())
    Br_spline = CubicSplineInterpolation((R_coords, Z_coords), Br_data; extrapolation_bc=Flat())
    Bz_spline = CubicSplineInterpolation((R_coords, Z_coords), Bz_data; extrapolation_bc=Flat())
    Bϕ_spline = CubicSplineInterpolation((R_coords, Z_coords), Bϕ_data; extrapolation_bc=Flat())

    return Plasma(
        R_coords,
        Z_coords,
        ne_spline,
        Br_spline,
        Bz_spline,
        Bϕ_spline)
end

# julia AD functions
dr_dτ(p, r, ϕ, z, kr, nϕ, kz, ω) = derivative(kr -> dispersion_relation(p, r, ϕ, z, kr, nϕ, kz, ω), kr)
dϕ_dτ(p, r, ϕ, z, kr, nϕ, kz, ω) = derivative(nϕ -> dispersion_relation(p, r, ϕ, z, kr, nϕ, kz, ω), nϕ)
dz_dτ(p, r, ϕ, z, kr, nϕ, kz, ω) = derivative(kz -> dispersion_relation(p, r, ϕ, z, kr, nϕ, kz, ω), kz)
dkr_dτ(p, r, ϕ, z, kr, nϕ, kz, ω) = derivative(r -> -dispersion_relation(p, r, ϕ, z, kr, nϕ, kz, ω), r)
dnϕ_dτ(p, r, ϕ, z, kr, nϕ, kz, ω) = derivative(ϕ -> -dispersion_relation(p, r, ϕ, z, kr, nϕ, kz, ω), ϕ)
dkz_dτ(p, r, ϕ, z, kr, nϕ, kz, ω) = derivative(z -> -dispersion_relation(p, r, ϕ, z, kr, nϕ, kz, ω), z)
dt_dτ(p, r, ϕ, z, kr, nϕ, kz, ω) = derivative(ω -> -dispersion_relation(p, r, ϕ, z, kr, nϕ, kz, ω), ω)

# Define the variables
@parameters τ ω p
@variables r(τ) ϕ(τ) z(τ) kr(τ) nϕ(τ) kz(τ) t(τ)
@register dr_dτ(p, r, ϕ, z, kr, nϕ, kz, ω)
@register dϕ_dτ(p, r, ϕ, z, kr, nϕ, kz, ω)
@register dz_dτ(p, r, ϕ, z, kr, nϕ, kz, ω)
@register dkr_dτ(p, r, ϕ, z, kr, nϕ, kz, ω)
@register dnϕ_dτ(p, r, ϕ, z, kr, nϕ, kz, ω)
@register dkz_dτ(p, r, ϕ, z, kr, nϕ, kz, ω)
@register dt_dτ(p, r, ϕ, z, kr, nϕ, kz, ω)

# System of ODEs using the AD functions
eqs = [
    Differential(τ)(r) ~ dr_dτ(p, r, ϕ, z, kr, nϕ, kz, ω),
    Differential(τ)(ϕ) ~ dϕ_dτ(p, r, ϕ, z, kr, nϕ, kz, ω),
    Differential(τ)(z) ~ dz_dτ(p, r, ϕ, z, kr, nϕ, kz, ω),
    Differential(τ)(kr) ~ dkr_dτ(p, r, ϕ, z, kr, nϕ, kz, ω),
    Differential(τ)(nϕ) ~ dnϕ_dτ(p, r, ϕ, z, kr, nϕ, kz, ω),
    Differential(τ)(kz) ~ dkz_dτ(p, r, ϕ, z, kr, nϕ, kz, ω),
    Differential(τ)(t) ~ dt_dτ(p, r, ϕ, z, kr, nϕ, kz, ω)
]

# Create an ODESystem object
@named ray_model = ModelingToolkit.ODESystem(eqs)

function solve(plasma::Plasma, r0::T, ϕ0::T, z0::T, θ_injection::T, freq::T, tmax::Float64) where {T<:Real}
    ω0 = 2π * freq

    # Define the parameters dictionary, if you have any parameters
    params = [ω => ω0, p => plasma]

    # figure out kr0 and kz0 based on injection angle
    res = Optim.optimize(kz -> dispersion_relation(plasma, r0, ϕ0, z0, 0.0, 0.0, kz, ω0), 0.0, 1E6, Optim.GoldenSection(), rel_tol=1E-3)
    k = res.minimizer
    kr0 = k * cos(θ_injection)
    kz0 = k * sin(θ_injection)
    nϕ0 = 0.0

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

function dispersion_relation(plasma::Plasma, r, ϕ, z, kr, nϕ, kz, ω)
    c = IMAS.constants.c

    kpar2, kperp2 = k_par2_perp2(plasma, r, ϕ, z, kr, nϕ, kz)

    S = S_stix(plasma, r, z, ω)
    D = D_stix(plasma, r, z, ω)
    P = P_stix(plasma, z, z, ω)

    term1 = S * (kperp2^2 * c^4) / ω^4
    term2 = ((S + P) * (S - (kpar2 * c^2) / ω^2) - D^2) * (kperp2 * c^2) / ω^2
    term3 = P * ((S - (kpar2 * c^2) / ω^2)^2 - D^2)

    return term1 - term2 + term3
end

function k_par2_perp2(plasma::Plasma, r, ϕ, z, kr, nϕ, kz)
    Br = plasma.Br_spline(r, z)
    Bϕ = plasma.Bϕ_spline(r, z)
    Bz = plasma.Bz_spline(r, z)
    B = sqrt(Br^2 + Bϕ^2 + Bz^2)
    kϕ = nϕ / r
    kpar2 = ((kr * Br + kz * Bz + kϕ * Bϕ) / B)^2.0

    k2 = kr^2 + kz^2 + kϕ^2
    kperp2 = k2 - kpar2

    return (kpar2=kpar2, kperp2=kperp2)
end

function S_stix(plasma::Plasma, r, z, ω)
    ω_pe = IMAS.ω_pe(plasma.ne_spline(r, z))
    ω_ce = IMAS.ω_ce(B_spline(plasma, r, z))
    return 1.0 - ω_pe^2 / (ω^2 - ω_ce^2)
end

function D_stix(plasma::Plasma, r, z, ω)
    ω_pe = IMAS.ω_pe(plasma.ne_spline(r, z))
    ω_ce = IMAS.ω_ce(B_spline(plasma, r, z))
    return ω_pe^2 * ω_ce / (ω * (ω^2 - ω_ce^2))
end

function P_stix(plasma::Plasma, r, z, ω)
    ω_pe = IMAS.ω_pe(plasma.ne_spline(r, z))
    return 1.0 - ω_pe^2 / ω^2
end

end # module TorJ
