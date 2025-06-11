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

ecrad_ref = open(dirname(@__DIR__) * "/../samples/ECRad_params.json", "r") do io
    JSON.parse(IOBuffer(read(io)))
end
ecrad_ref["s"] = Vector{Float64}(ecrad_ref["s"])
ecrad_ref["R"] = Vector{Float64}(ecrad_ref["R"])
ecrad_ref["z"] = Vector{Float64}(ecrad_ref["z"])

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
@time u = TorJ.make_ray(plasma, x0, y0, z0, β, α, freq, 1, 0.4);

x = hcat([tu[1:3] for tu in u]...)
R = hypot.(x[1,:], x[2,:])
R0 = min(ecrad_ref["R"][1], R[1])
s_torj = zeros(Float64, size(R))
s_torj[2:end] = hypot.(R[2:end] .- R[1:(end - 1)], x[3, 2:end] .- x[3, 1:(end - 1)])
s_torj .= cumsum(s_torj)
i_ref = 1
for i in eachindex(ecrad_ref["R"])
    global i_ref = i
    if abs(ecrad_ref["R"][i] - R0) < abs(ecrad_ref["R"][i+1] - R0) break end
end
s_offset_ecrad = ecrad_ref["s"][i_ref]
i_ref = 1
for i in eachindex(ecrad_ref["R"])
    global i_ref = i
    if abs(R[i] - R0) < abs(R[i+1] - R0) break end
end
s_offset_torj = s_torj[i_ref]
s_torj .-= s_offset_torj
ecrad_r_itp = IMAS.interp1d(ecrad_ref["s"] .- s_offset_ecrad, ecrad_ref["R"], :cubic)
ecrad_z_itp = IMAS.interp1d(ecrad_ref["s"] .- s_offset_ecrad, ecrad_ref["z"], :cubic)
for i in eachindex(R)
    psi_norm = TorJ.evaluate(plasma.psi_norm_spline, x[:,i])
    if psi_norm >= maximum(ecrad_ref["rhop"])^2 continue end
    dist = hypot.(ecrad_r_itp(s_torj[i]) .- R[i], ecrad_z_itp(s_torj[i]) .- x[3,i])
    @test dist < 2.e-3
end