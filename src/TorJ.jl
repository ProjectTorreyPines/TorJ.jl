module TorJ

import Interpolations: CubicSplineInterpolation, Line, Extrapolation, gradient
import LinearAlgebra
import IMAS
using DifferentialEquations
using NLsolve
using Roots
import Enzyme: autodiff, Duplicated, Const, Reverse

include("plasma.jl")

include("dispersion.jl")

include("solve.jl")


end # module TorJ
