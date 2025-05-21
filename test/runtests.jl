import Pkg
# Pkg.activate(@__DIR__)

using Test
import IMAS
import TorJ


dd = IMAS.json2imas(dirname(@__DIR__) *"/samples/D3D_170325_trimmed.json"; error_on_missing_coordinates=false)
dd.global_time = 2.0
eq_slice = dd.equilibrium.time_slice[]
R_grid = eq_slice.profiles_2d[1].grid.dim1
z_grid = eq_slice.profiles_2d[1].grid.dim2

profiles_1d = dd.core_profiles.profiles_1d[]


plasma = TorJ.Plasma(R_grid, z_grid, eq_slice.profiles_2d[1].psi, 
                     profiles_1d.grid.psi, profiles_1d.electrons.density,
                     eq_slice.profiles_2d[1].b_field_r, eq_slice.profiles_2d[1].b_field_z,
                     eq_slice.profiles_2d[1].b_field_tor);

freq = 1.1E11
r0 = 4.30
z0 = 0.3
θ0 = 0.0
ϕ0 = 0.0
t0 = 0.0
@time sol = TorJ.make_ray(plasma, r0, 0.0, z0, ϕ0, θ0, freq, -1, 1E2);
