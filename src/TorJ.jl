module TorJ

import Optim
import Interpolations: CubicSplineInterpolation, Flat, Extrapolation, gradient
import ForwardDiff
import LinearAlgebra
import IMAS
using DifferentialEquations
using NLsolve
using Roots

include("plasma.jl")

include("dispersion.jl")

include("solve.jl")

end # module TorJ
