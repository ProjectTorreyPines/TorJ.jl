include("setup.jl")
#Higher frequency for full beam test to avoid cut-off
# Check Dagger setup
println("Checking Dagger setup...")
TorJ.abs_Al_init(31)
dP_dV_psi = Vector(LinRange(0.0, 1.0, 1000))
@testset "make_beam test" begin
    @time arclengths, trajectories, ray_powers, dP_dV, ray_weights, absorbed_power_fraction = TorJ.make_beam(plasma_low_density, R0, phi0, z0, steering_angle_tor, 
                                        steering_angle_pol, spot_size, 
                                        inverse_curvature_radius, f_abs_test, 1, 1.0, dP_dV_psi);
    # Compare total power absorbed against TORBEAM
    # This case has nearly total absorption but not quite making it a decent benchmark
    println("Comparing TorJ P_abs with Torbeam: $absorbed_power_fraction, $(tb_ref["P_abs"]/1.e6)")
    @test isapprox(absorbed_power_fraction, tb_ref["P_abs"]/1.e6; atol=0.001, rtol=1.e-3)

    absorbed_power = 0.0
    for i_ray in eachindex(ray_weights)
        absorbed_power -= ray_powers[i_ray][end]*ray_weights[i_ray]
    end
    absorbed_power += 1.0
    # Check if power deposition profile is consistent with results from ray tracing
    println("Checking power deposition and ray absorbed power consistency: $absorbed_power_fraction, $absorbed_power")
    @test isapprox(absorbed_power_fraction, absorbed_power; atol=0.001, rtol=1.e-3)
end
