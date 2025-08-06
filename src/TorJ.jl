module TorJ

import Optim
import Interpolations: cubic_spline_interpolation, Line, Extrapolation, gradient
import ForwardDiff
import LinearAlgebra
import IMAS
using DifferentialEquations
using NLsolve
using Roots
using SpecialFunctions: besselj
using FastGaussQuadrature
using Dagger
import Dierckx

# Module-level constants (equivalent to Fortran parameters)


include("constants.jl")

include("launch.jl")

include("plasma.jl")

include("dispersion.jl")

include("solve.jl")

include("absorption.jl")

end # module TorJ
