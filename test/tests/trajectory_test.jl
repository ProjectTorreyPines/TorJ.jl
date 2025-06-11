using Test
import IMAS
import TorJ
using Plots
using JSON

dd = IMAS.json2imas(dirname(@__DIR__) *"/../samples/sample_L_mode_peaked.json"; error_on_missing_coordinates=false)
dd.global_time = 2.0
eq_slice = dd.equilibrium.time_slice[]
R_grid = eq_slice.profiles_2d[].grid.dim1
z_grid = eq_slice.profiles_2d[].grid.dim2

profiles_1d = dd.core_profiles.profiles_1d[]
freq = 85.5E9
# (eq_slice.profiles_2d[1].psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis), 
# (profiles_1d.grid.psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis), 
psi_norm_ne = (profiles_1d.grid.psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis) 
plasma = TorJ.Plasma(R_grid, z_grid, 
                     (eq_slice.profiles_2d[1].psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis),
                      psi_norm_ne,
                      profiles_1d.electrons.density,
                      eq_slice.profiles_2d[1].b_field_r, eq_slice.profiles_2d[1].b_field_z,
                      eq_slice.profiles_2d[1].b_field_tor);

data = Dict("rho_tor_norm" => profiles_1d.grid.rho_tor_norm, "psi_norm" => psi_norm_ne, "n_e" => profiles_1d.electrons.density)
json_str = JSON.json(data)
open(dirname(@__DIR__) *"/tests/profile.json", "w") do io
    write(io, json_str)
end

ecrad_ref = open(dirname(@__DIR__) * "/../samples/ECRad_params.json", "r") do io
    JSON.parse(IOBuffer(read(io)))
end

for i in 1:length(ecrad_ref["x"])
    x = vec([ecrad_ref["x"][i]; ecrad_ref["y"][i]; ecrad_ref["z"][i]])
    N = vec([ecrad_ref["Nx"][i]; ecrad_ref["Ny"][i]; ecrad_ref["Nz"][i]])
    B_test = TorJ.B_spline(plasma,x)
    B_diff = B_test .- vec([ecrad_ref["Bx"][i];ecrad_ref["By"][i];ecrad_ref["Bz"][i]])
    @test maximum(abs.(B_diff)) < 1.e-6
    ne = TorJ.n_e(plasma, x)
    ne_rel_diff = (ne - ecrad_ref["ne"][i]) / (2.0 * (ne + ecrad_ref["ne"][i]))
    @test abs(ne_rel_diff) < 1.e-2
    X, Y = TorJ.eval_plasma(plasma, x, N, 2.0 * pi * freq)
    Y_diff = abs(Y - ecrad_ref["Y"][i])
    @test Y_diff < 1.e-6
end