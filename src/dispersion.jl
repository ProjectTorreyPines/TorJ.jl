"""
    dispersion_relation(plasma::Plasma, r, ϕ, z, kr, nϕ, kz, ω)

Cold plasma dispersion relation
"""

function eval_plasma(plasma::Plasma, x::AbstractVector{<:Real}, N::AbstractVector{<:Real}, omega::Real)
    B = B_spline(plasma, x)
    B_abs = LinearAlgebra.norm(B)
    b = B ./ B_abs
    N_par = LinearAlgebra.dot(N, b) 
    X = n_e(plasma, x) * constants.e^2/(constants.ϵ_0 * constants.m_e * omega^2)
    Y = constants.e * B_abs / (constants.m_e * omega)
    return X, Y, N_par, b
end

function wrap_eval_plasma(plasma::Plasma, x::AbstractVector{<:Real}, N::AbstractVector{<:Real}, omega::Real, index::Integer)
    return eval_plasma(plasma, x, N, omega)[index]
end

function Δ( X::Real, Y::Real, N_par::Real)
    return (1.0 - N_par^2)^2 + 4.0 * N_par^2 * (1.0 - X) / Y^2
end

function sqrtΔ( X::Real, Y::Real, N_par::Real)
    return sqrt(Δ(X, Y, N_par))
end

function refractive_index_sq( X::Real, Y::Real, N_par::Real, mode::Integer)
    Ns_sq = 1.e0 - X + (1.0 + Real(mode) * sqrtΔ(X, Y, N_par) + N_par^2)/(2.0 * (-1.0 + X + Y^2)) * X * Y^2
    return Ns_sq
end

function dispersion_relation(x::AbstractVector{<:Real},  N::AbstractVector{<:Real}, plasma::Plasma, omega:: Real, mode::Integer)
    N_abs = LinearAlgebra.norm(N)
    X, Y, N_par, b = eval_plasma(plasma, x, N, omega)
    Ns_sq = refractive_index_sq(X, Y, N_par, mode)
    return N_abs^2 - Ns_sq
end

function dispersion_relation_enzyme(u::AbstractVector{<:Real},  plasma::Plasma, omega:: Real, mode::Integer)
    x = @view u[1:3]
    N = @view u[4:6]
    return dispersion_relation(x, N, plasma, omega, mode)
end

function dΛ_dN_ana(N:: AbstractVector{<:Real}, X::Real, Y::Real, N_par::Real, b::AbstractVector{<:Real}, mode::Integer)
    dNs_sq_dN_par = (N_par*X*Y^2*(1.0 + (Real(mode)*(2.0 - 2*X + (-1.0 + N_par^2)*Y^2))/ (Y^2*sqrtΔ(X, Y, N_par))))/(-1.e0 + X + Y^2)
    return 2.e0 * N .- b .* dNs_sq_dN_par
end


