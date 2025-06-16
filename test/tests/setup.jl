using Test
import IMAS
import TorJ
import Artifacts
import Pkg
using JSON

# Only run preprocessing if it hasn't been done already
if !@isdefined(TEST_DATA_LOADED)

    Pkg.ensure_artifact_installed("data", "Artifacts.toml")
    artifact_path = Artifacts.artifact"data"

    ecrad_ref = open(artifact_path * "/data/ECRad_params.json", "r") do io
        JSON.parse(IOBuffer(read(io)))
    end

    dd = IMAS.json2imas(artifact_path *"/data/sample_L_mode_peaked.json"; error_on_missing_coordinates=false)
    dd.global_time = 2.0

    eq_slice = dd.equilibrium.time_slice[]
    R_grid = eq_slice.profiles_2d[].grid.dim1
    z_grid = eq_slice.profiles_2d[].grid.dim2

    profiles_1d = dd.core_profiles.profiles_1d[]


    ecrad_ref["s"] = Vector{Float64}(ecrad_ref["s"])
    ecrad_ref["R"] = Vector{Float64}(ecrad_ref["R"])
    ecrad_ref["z"] = Vector{Float64}(ecrad_ref["z"])

    plasma = TorJ.Plasma(R_grid, z_grid, 
                        (eq_slice.profiles_2d[1].psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis),
                        (profiles_1d.grid.psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis),
                        profiles_1d.electrons.density,
                        profiles_1d.electrons.temperature,
                        eq_slice.profiles_2d[1].b_field_r, eq_slice.profiles_2d[1].b_field_z,
                        eq_slice.profiles_2d[1].b_field_tor);

    freq = 85.5E9
    x0 = 2.5
    y0 = 0.0
    z0 = 0.4
    α = deg2rad(30.0) # Convert from TORBEAM convetion to IMAS (they are the same for phi_tor ==0)
    β = 0.0
    global TEST_DATA_LOADED = true
end