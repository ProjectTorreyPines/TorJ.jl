using Test
import TorJ
import TorJ: IMAS
import Artifacts
import Pkg
import JSON
if !@isdefined(FORCE_RELOAD_TEST_DATA)
    FORCE_RELOAD_TEST_DATA = false
end
# Only run preprocessing if it hasn't been done already
if !@isdefined(TEST_DATA_LOADED) || FORCE_RELOAD_TEST_DATA
    Pkg.ensure_artifact_installed("data", joinpath(@__DIR__, "..", "Artifacts.toml"))
    artifact_path = Artifacts.artifact"data"

    ecrad_ref = open(artifact_path * "/data/ECRad_params.json", "r") do io
        JSON.parse(IOBuffer(read(io)))
    end

    ecrad_ref_abs = open(artifact_path * "/data/ECRad_params_2.json", "r") do io
        JSON.parse(IOBuffer(read(io)))
    end

    tb_ref = open(artifact_path * "/data/tb_results.json", "r") do io
        JSON.parse(IOBuffer(read(io)))
    end

    toray_ref = open(artifact_path * "/data/toray_results.json", "r") do io
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
                        eq_slice.profiles_2d[1].b_field_tor,
                        eq_slice.profiles_1d.psi, eq_slice.profiles_1d.volume);
    # For comparing against TORBEAM we want less dispersion
    plasma_low_density = TorJ.Plasma(R_grid, z_grid,
                        (eq_slice.profiles_2d[1].psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis),
                        (profiles_1d.grid.psi .- eq_slice.global_quantities.psi_axis)./(eq_slice.global_quantities.psi_boundary - eq_slice.global_quantities.psi_axis),
                        profiles_1d.electrons.density * 0.3,
                        profiles_1d.electrons.temperature,
                        eq_slice.profiles_2d[1].b_field_r, eq_slice.profiles_2d[1].b_field_z,
                        eq_slice.profiles_2d[1].b_field_tor,
                        eq_slice.profiles_1d.psi, eq_slice.profiles_1d.volume);

    f = 85.5E9
    f_abs_test = 92.5E9
    R0 = 2.5
    phi0 = 0.0
    x0 = R0 * cos(phi0)
    y0 = R0 * sin(phi0)
    z0 = 0.4
    spot_size = 0.0174  # beam width parameter
    inverse_curvature_radius = 1.0/3.99
    steering_angle_pol = deg2rad(30.0) # Convert from TORBEAM convetion to IMAS (they are the same for phi_tor ==0)
    steering_angle_tor = 0.0

    # Define psi grid for dP_dV calculation
    psi_dP_dV = Vector(LinRange(0.0, 1.0, 1000))

    # Set up the points for the resonance ellipse integration in the absorption coefficient
    TorJ.abs_Al_init(24)
    # Set this global so we don't repeat this when we `include` it again
    global TEST_DATA_LOADED = true
end