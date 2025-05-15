module TorJ

import Optim
import Interpolations: CubicSplineInterpolation, Flat, Extrapolation
import ForwardDiff: derivative
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: Symbolics
using IMAS
using NLsolve


include("plasma.jl")

include("dispersion.jl")

include("solve.jl")

end # module TorJ
