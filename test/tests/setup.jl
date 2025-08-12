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

    plasma = TorJ.Plasma(dd);
    # For comparing against TORBEAM we want less dispersion
    plasma_low_density = TorJ.Plasma(dd; ne_scale=0.3);

    ecrad_ref["s"] = Vector{Float64}(ecrad_ref["s"])
    ecrad_ref["R"] = Vector{Float64}(ecrad_ref["R"])
    ecrad_ref["z"] = Vector{Float64}(ecrad_ref["z"])


    f = 85.5E9
    
    R0 = 2.5
    phi0 = 0.0
    x0 = R0 * cos(phi0)
    y0 = R0 * sin(phi0)
    z0 = 0.4
    spot_size = 0.0174  # beam width parameter
    inverse_curvature_radius = 1.0/3.99
    steering_angle_pol = deg2rad(30.0) # Convert from TORBEAM convetion to IMAS (they are the same for phi_tor ==0)
    steering_angle_tor = 0.0

    f_abs_test = 92.5E9
    # Define psi grid for dP_dV calculation
    psi_dP_dV = Vector(LinRange(0.0, 1.0, 1000))

    # Set up the points for the resonance ellipse integration in the absorption coefficient
    TorJ.abs_Al_init(24)

    # Fill in ec_launchers and pulse_schedule
    IMAS.resize!(dd.ec_launchers.beam, 2)
    dd.ec_launchers.beam[1].time = ones(Float64, 1) .* 2.0
    dd.ec_launchers.beam[1].frequency.data = ones(Float64, 1) .* f
    dd.ec_launchers.beam[1].frequency.time = ones(Float64, 1) .* 2.0
    dd.ec_launchers.beam[1].power_launched.data = ones(Float64, 1) .* 1.e6
    dd.ec_launchers.beam[1].power_launched.time = ones(Float64, 1) .* 2.0
    dd.ec_launchers.beam[1].mode = -1

    dd.ec_launchers.beam[1].launching_position.r =  ones(Float64, 1) .* R0
    dd.ec_launchers.beam[1].launching_position.phi = ones(Float64, 1) .* phi0
    dd.ec_launchers.beam[1].launching_position.z = ones(Float64, 1) .* z0

    dd.ec_launchers.beam[1].steering_angle_pol = ones(Float64, 1) .* steering_angle_pol
    dd.ec_launchers.beam[1].steering_angle_tor = ones(Float64, 1) .* steering_angle_tor

    dd.ec_launchers.beam[1].spot.size = ones(Float64, (2, 1)) .* spot_size
    dd.ec_launchers.beam[1].phase.curvature = ones(Float64, (2, 1)) .* inverse_curvature_radius

    dd.ec_launchers.beam[2] = deepcopy(dd.ec_launchers.beam[1])
    dd.ec_launchers.beam[2].frequency.data[:] .= f_abs_test

    IMAS.resize!(dd.pulse_schedule.ec.beam, 2)
    dd.pulse_schedule.ec.time = ones(Float64, 1) .* 2.0
    dd.pulse_schedule.ec.beam[1].power_launched.reference = ones(Float64, 1) .* 1.e6
    dd.pulse_schedule.ec.beam[2].power_launched.reference = ones(Float64, 1) .* 1.e6

    # Set this global so we don't repeat this when we `include` it again
    global TEST_DATA_LOADED = true
end