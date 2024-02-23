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
@named ray_model = ModelingToolkit.ODESystem(eqs)#; discrete_events = [ (r > 9.0) => integrator -> terminate!(integrator)])