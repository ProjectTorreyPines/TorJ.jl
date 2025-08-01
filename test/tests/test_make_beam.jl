include("setup.jl")
@testset "Raytrace test" begin
    @time trajectories = TorJ.make_beam(plasma, R0, phi0, z0, β, α, 
            spot_size, inverse_curvature_radius, freq, 1, 0.4);
end