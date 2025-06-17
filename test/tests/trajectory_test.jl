
include("setup.jl")
@testset "Trajectory parameters test" begin

    for i in 1:length(ecrad_ref["x"])
        x = [ecrad_ref["x"][i] ecrad_ref["y"][i]; ecrad_ref["z"][i]]
        N = [ecrad_ref["Nx"][i], ecrad_ref["Ny"][i], ecrad_ref["Nz"][i]]
        B_test = TorJ.B_spline(plasma,x)
        B_diff = B_test .- [ecrad_ref["Bx"][i], ecrad_ref["By"][i], ecrad_ref["Bz"][i]]
        @test maximum(abs.(B_diff)) < 1.e-6
        ne = TorJ.n_e(plasma, x)
        ne_rel_diff = (ne - ecrad_ref["ne"][i]) / (2.0 * (ne + ecrad_ref["ne"][i]))
        @test abs(ne_rel_diff) < 1.e-2
        X, Y = TorJ.eval_plasma(plasma, x, N, 2.0 * pi * freq)
        Y_diff = abs(Y - ecrad_ref["Y"][i])
        @test Y_diff < 1.e-6
    end
end