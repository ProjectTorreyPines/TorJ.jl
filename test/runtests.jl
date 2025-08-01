if !isempty(ARGS)
    # Run only specific tests
    for test_name in ARGS
        include("tests/$(test_name)_test.jl")
    end
else
    include("tests/raytrace_test.jl")
    include("tests/trajectory_test.jl")
    include("tests/absorption_test.jl")
    include("tests/test_launch_weights.jl")
end


