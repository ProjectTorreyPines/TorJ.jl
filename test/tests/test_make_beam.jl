include("setup.jl")

# Check Dagger setup
println("Checking Dagger setup...")
TorJ.check_dagger_setup()

@testset "Raytrace test" begin
    @time arclengths, trajectories, ray_powers, ray_weights = TorJ.make_beam(plasma, R0, phi0, z0, steering_angle_tor, 
                                        steering_angle_pol, spot_size, 
                                        inverse_curvature_radius, freq, 1, 0.4);
end