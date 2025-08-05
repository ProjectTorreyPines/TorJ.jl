include("setup.jl")
#Higher frequency for full beam test to avoid cut-off
# Check Dagger setup
println("Checking Dagger setup...")
TorJ.abs_Al_init(31)
dP_dV_psi = LinRange(0.0, 1.0, 1000)
@testset "process_launcher test" begin
    @time arclengths, trajectories, ray_powers, dP_dV, ray_weights, absorbed_power_fraction = TorJ.process_launcher(plasma, R0, phi0, z0, steering_angle_tor, 
                                        steering_angle_pol, spot_size, 
                                        inverse_curvature_radius, f_abs_test, 1, 1.0, dP_dV_psi);
    # Compare total power absorbed against TORBEAM
    # This case has nearly total absorption but not quite making it a decent benchmark
    @test isapprox(absorbed_power_fraction, tb_ref["P_abs"]; atol=0.001, rtol=1.e-3)
end