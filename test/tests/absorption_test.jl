include("setup.jl")

@testset "Absorption test" begin
    omega = 2.0 * pi * freq
    TorJ.set_extv!()
    α = zeros(Float64, length(ecrad_ref["x"]))
    Nr = ecrad_ref["Nc"][1]
    for i in 1:length(ecrad_ref["x"])
        α[i] = TorJ.abs_Albajar_fast(omega, ecrad_ref["X"][i], ecrad_ref["Y"][i], 
                                     ecrad_ref["theta"][i], ecrad_ref["Te"][i], 1)
    end
    if !all(α .≈ ecrad_ref["ab"])
        mask = .!(α .≈ ecrad_ref["ab"])
        for i in 1:length(α[mask])
            println(α[mask][i], " ", ecrad_ref["ab"][mask][i])
        end
    end
    @test all(α .≈ ecrad_ref["ab"])
end