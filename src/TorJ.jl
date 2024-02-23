module TorJ

import Optim
import Interpolations: CubicSplineInterpolation, Flat, AbstractInterpolation
import ForwardDiff: derivative
using ModelingToolkit, OrdinaryDiffEq

include("plasma.jl")

include("stix.jl")

include("dispersion.jl")

include("ode.jl")

include("solve.jl")

end # module TorJ
