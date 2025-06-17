include("setup.jl")

@testset "Absorption test" begin
    TorJ.set_extv!()
    for i in 1:length(ecrad_ref["x"])
        α = TorJ.α(ecrad_ref["X"][i], ecrad_ref["Y"][i], ecrad_ref["Nc"][i], 
                              ecrad_ref["theta"][i], ecrad_ref["Te"][i], ecrad_ref["v_g_perp"][i], 1)
        @test α == ecrad_ref["ab_second"][i]
    end
end