# julia AD functions
dx_ds(p, x, y, z, Nx, Ny, Nz, ω, mode) = derivative(kx -> dispersion_relation(p, x, y, z, Nx, Ny, Nz, ω, mode), kx)
dy_ds(p, x, y, z, Nx, Ny, Nz, ω, mode) = derivative(ky -> dispersion_relation(p, x, y, z, Nx, Ny, Nz, ω, mode), ky)
dz_ds(p, x, y, z, Nx, Ny, Nz, ω, mode) = derivative(kz -> dispersion_relation(p, x, y, z, Nx, Ny, Nz, ω, mode), kz)
dNx_ds(p, x, y, z, Nx, Ny, Nz, ω, mode) = derivative(x -> -dispersion_relation(p, x, y, z, Nx, Ny, Nz, ω, mode), x)
dNy_ds(p, x, y, z, Nx, Ny, Nz, ω, mode) = derivative(y -> -dispersion_relation(p, x, y, z, Nx, Ny, Nz, ω, mode), y)
dNz_ds(p, x, y, z, Nx, Ny, Nz, ω, mode) = derivative(z -> -dispersion_relation(p, x, y, z, Nx, Ny, Nz, ω, mode), z)

# Normalizing function
ds_dτ(p, x, y, z, Nx, Ny, Nz, ω, mode) = sqrt(dx_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)^2 +
                                              dy_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)^2 +
                                              dz_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)^2)

# Define the variables
@parameters s ω p
@variables x(s) y(s) z(s) Nx(s) Ny(s) Nz(s)
@register dx_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)
@register dy_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)
@register dz_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)
@register dNx_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)
@register dNy_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)
@register dNz_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)
@register ds_dτ(p, x, y, z, Nx, Ny, Nz, ω, mode)

# System of ODEs using the AD functions
eqs = [
    Differential(s)(x) ~ dx_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)/ds_dτ(p, x, y, z, Nx, Ny, Nz, ω, mode),
    Differential(s)(y) ~ dy_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)/ds_dτ(p, x, y, z, Nx, Ny, Nz, ω, mode),
    Differential(s)(z) ~ dz_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)/ds_dτ(p, x, y, z, Nx, Ny, Nz, ω, mode),
    Differential(s)(Nx) ~ dkr_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)/ds_dτ(p, x, y, z, Nx, Ny, Nz, ω, mode),
    Differential(s)(Ny) ~ dnϕ_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)/ds_dτ(p, x, y, z, Nx, Ny, Nz, ω, mode),
    Differential(s)(Nz) ~ dkz_ds(p, x, y, z, Nx, Ny, Nz, ω, mode)/ds_dτ(p, x, y, z, Nx, Ny, Nz, ω, mode),
]
# Create an ODESystem object
@named ray_model = ModelingToolkit.ODESystem(eqs)#; discrete_events = [ (r > 9.0) => integrator -> terminate!(integrator)])