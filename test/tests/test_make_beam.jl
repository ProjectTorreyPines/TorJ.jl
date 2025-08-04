include("setup.jl")
#Higher frequency for full beam test to avoid cut-off
# Check Dagger setup
println("Checking Dagger setup...")
TorJ.abs_Al_init(31)
TorJ.check_dagger_setup()

@testset "make_beam test" begin
    @time arclengths, trajectories, ray_powers, ray_weights = TorJ.make_beam(plasma, R0, phi0, z0, steering_angle_tor, 
                                        steering_angle_pol, spot_size, 
                                        inverse_curvature_radius, f_abs_test, 1, 0.4);
end