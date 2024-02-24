import Pkg
Pkg.activate(@__DIR__)

using Test

import TorJ
import BSON

BSON.@load joinpath(@__DIR__, "TorJsample_data.bson") R_coords Z_coords ne_data Br_data Bz_data Bϕ_data

plasma = TorJ.Plasma(R_coords, Z_coords, ne_data, Br_data, Bz_data, Bϕ_data);

freq = 4.5E9
r0 = 8.30
ϕ0 = 0.0
z0 = 0.8
θ0 = -1.0 * π
nϕ0 = 0.0
t0 = 0.0
ω0 = 2π * freq

@time sol = TorJ.solve(plasma, r0, ϕ0, z0, nϕ0, θ0, freq, 1E2);
