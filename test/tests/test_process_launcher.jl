include("setup.jl")
#Higher frequency for full beam test to avoid cut-off
# Check Dagger setup
println("Checking Dagger setup...")
TorJ.abs_Al_init(31)
TorJ.check_dagger_setup()
dP_dV_psi = LinRange(0.0, 1.0, 1000)
@testset "process_launcher test" begin
    @time arclengths, trajectories, ray_powers, ray_weights, dP_dV = TorJ.process_launcher(plasma, R0, phi0, z0, steering_angle_tor, 
                                        steering_angle_pol, spot_size, 
                                        inverse_curvature_radius, f_abs_test, 1, 1.0, dP_dV_psi);
end