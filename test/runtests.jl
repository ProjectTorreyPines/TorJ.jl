if !isempty(ARGS)
    # Run only specific tests
    for test_name in ARGS
        include("tests/$(test_name)_test.jl")
    end
else
    include("tests/test_absorption.jl")
    include("tests/test_make_ray.jl")
    include("tests/test_trajectory.jl")
    include("tests/test_absorption.jl")
    include("tests/test_launch_weights.jl")
    include("tests/test_make_beam.jl")
    include("tests/test_process_ids.jl")
end


