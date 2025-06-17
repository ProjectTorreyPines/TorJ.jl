include("setup.jl")

@testset "Absorption test" begin
    TorJ.set_extv!()
    for i in 1:length(ecrad_ref["x"])
        @test TorJ.α(ecrad_ref["X"][i], ecrad_ref["Y"][i], ecrad_ref["Nc"][i], 
                              ecrad_ref["theta"][i], ecrad_ref["Te"][i], 1) == ecrad_ref["ab_second"][i]
    end
end