include("setup.jl")

@testset "Absorption test" begin
    omega = 2.0 * pi * freq_abs_test
    TorJ.abs_Al_init(24)
    α = zeros(Float64, length(ecrad_ref_abs["x"]))
    Nr = ecrad_ref_abs["Nc"][1]
    for i in 1:length(ecrad_ref_abs["x"])
        if i == 1
            ds = ecrad_ref_abs["s"][2] - ecrad_ref_abs["s"][1]
        else
            ds = ecrad_ref_abs["s"][i] - ecrad_ref_abs["s"][i-1]
        end
        α[i] = TorJ.abs_Albajar_fast(omega, ecrad_ref_abs["X"][i], ecrad_ref_abs["Y"][i], 
                                     ecrad_ref_abs["theta"][i], ecrad_ref_abs["Te"][i], 1, abs(ds))
    end
    if !all(isapprox(α, ecrad_ref_abs["ab"];atol=0.1, rtol=1.e-2))
        mask = .!(isapprox(α, ecrad_ref_abs["ab"];atol=0.1, rtol=1.e-2))
        for i in 1:length(α[mask])
            println(α[mask][i], " ", ecrad_ref_abs["ab"][mask][i])
        end
    end
    @test all(isapprox(α, ecrad_ref_abs["ab"];atol=0.1, rtol=1.e-2))
end