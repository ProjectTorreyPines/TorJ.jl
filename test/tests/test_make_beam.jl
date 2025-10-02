include("setup.jl")
import TorJ: Dierckx
#Higher frequency for full beam test to avoid cut-off

@testset "make_beam test" begin
    println("Using $(Threads.nthreads()) threads")
    @time arclengths, trajectories, ray_powers, dP_dV, absorbed_power_fraction, ray_weights = TorJ.make_beam(plasma_low_density, R0, phi0, z0, steering_angle_tor, 
                                        steering_angle_pol, spot_size, 
                                        inverse_curvature_radius, f_abs_test, 1, 1.0, psi_dP_dV);
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
    # Check that the power deposition profile yields the absorbed power when integrated over the plasma volume
    eqt = dd.equilibrium.time_slice[]
    eqglobs = eqt.global_quantities
    eqt1d = eqt.profiles_1d
    _norm = psi -> (psi .- eqglobs.psi_axis) ./ (eqglobs.psi_boundary - eqglobs.psi_axis)
    volume_psi_spl = Dierckx.Spline1D(_norm(eqt1d.psi), eqt1d.volume)
    _norm(eqt1d.psi)
    # This is equidistant
    dpsi = psi_dP_dV[2] - psi_dP_dV[1]
    P_test = 0.0
    for (i, psi_norm) in enumerate(psi_dP_dV)
        P_test += Dierckx.derivative(volume_psi_spl, psi_norm, nu=1) * dP_dV[i] * dpsi
    end
    println("Checking manual volume integral against ray absorbed power for consistency: $absorbed_power_fraction, $P_test")
    @test isapprox(absorbed_power_fraction, P_test; atol=0.001, rtol=1.e-3)
end
