module TorJ

import Optim
import Interpolations: CubicSplineInterpolation, Line, Extrapolation, gradient
import ForwardDiff
import LinearAlgebra
import IMAS
using DifferentialEquations
using NLsolve
using Roots
using SpecialFunctions: sphericalbesselj


include("plasma.jl")

include("dispersion.jl")

include("solve.jl")

include("absorption.jl")
import .absorption: set_extv!, α

end # module TorJ
