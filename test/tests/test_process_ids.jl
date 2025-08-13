include("setup.jl")
# Simple functionality test.
@testset "Raytrace test" begin
    println("Using $(Threads.nthreads()) threads")
    TorJ.process_ids!(dd)
end