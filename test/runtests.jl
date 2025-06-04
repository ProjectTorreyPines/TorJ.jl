import Pkg
# Pkg.activate(@__DIR__)

using Test
import IMAS
import TorJ
using Plots
using JSON

dd = IMAS.json2imas(dirname(@__DIR__) *"/samples/sample_L_mode_peaked.json"; error_on_missing_coordinates=false)
dd.global_time = 2.0
eq_slice = dd.equilibrium.time_slice[]
R_grid = eq_slice.profiles_2d[].grid.dim1
z_grid = eq_slice.profiles_2d[].grid.dim2

profiles_1d = dd.core_profiles.profiles_1d[]


plasma = TorJ.Plasma(R_grid, z_grid, 
                    (eq_slice.profiles_2d[1].psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis), 
                    (profiles_1d.grid.psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis), 
                     profiles_1d.electrons.density,
                     eq_slice.profiles_2d[1].b_field_r, eq_slice.profiles_2d[1].b_field_z,
                     eq_slice.profiles_2d[1].b_field_tor);

freq = 85.5E9
x0 = 2.5
y0 = 0.0
z0 = 0.4
α = deg2rad(30.0) # Convert from TORBEAM convetion to IMAS (they are the same for phi_tor ==0)
β = 0.0
u = TorJ.make_ray(plasma, x0, y0, z0, β, α, freq, 1, 0.4);
x = hcat([tu[1:3] for tu in u]...)
R = hypot.(x[1,:], x[2,:])
p =plot(R, x[3,:]; xlim=(0.7, 2.3), ylim=(-1.0,1.0), aspect_ratio=:equal);
savefig(p, dirname(@__DIR__) *"/plots/first_ray.pdf")

data = Dict("R" => R, "z" => x[3,:])
json_str = JSON.json(data)

open(dirname(@__DIR__) *"/test/ray_trajectory.json", "w") do io
    write(io, json_str)
end

