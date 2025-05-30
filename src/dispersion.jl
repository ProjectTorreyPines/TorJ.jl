"""
    dispersion_relation(plasma::Plasma, r, ϕ, z, kr, nϕ, kz, ω)

Cold plasma dispersion relation
"""

const constants = (
    μ_0=1.25663706212e-6,
    c=2.99792458e8,
    ϵ_0=8.8541878128e-12,
    k_B=1.380649e-23,
    e=1.602176634e-19,
    m_e=9.1093837015e-31,
    m_p=1.67262192369e-27,
    m_n=1.67492749804e-27,
    atm=101325.0,
    m_u=1.6605390666e-27,
    avog=6.02214076e23
)


function eval_plasma(plasma::Plasma, x::AbstractVector{<:Real}, N::AbstractVector{<:Real}, omega::Real)
    r = hypot(x[1], x[2])
    phi = atan(x[2], x[1])
    B_cyl = B_spline(plasma, r, phi, x[3])
    B = zeros(typeof(phi), 3)
    B[1] = B_cyl[1] * cos(phi) - B_cyl[2] * sin(phi)
    B[2] = B_cyl[1] * sin(phi) + B_cyl[2] * cos(phi)
    B[3] = B_cyl[3]
    B_abs = LinearAlgebra.norm(B)
    b = B ./ B_abs
    N_par = LinearAlgebra.dot(N, b) 
    n_e = exp(plasma.log_ne_spline(r, x[3])) #
    X = n_e * constants.e^2/(constants.ϵ_0 * constants.m_e * omega^2)
    Y = constants.e* B_abs / (constants.m_e * omega)
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

function dispersion_relation(x::AbstractVector{<:Real}, N::AbstractVector{<:Real}, plasma::Plasma, omega:: Real, mode::Integer)
    N_abs = LinearAlgebra.norm(N)
    X, Y, N_par, b = eval_plasma(plasma, x, N, omega)
    Ns_sq = refractive_index_sq(X, Y, N_par, mode)
    return N_abs^2 - Ns_sq
end
function dΛ_dN_ana(N:: AbstractVector{<:Real}, X::Real, Y::Real, N_par::Real, b::AbstractVector{<:Real}, mode::Integer)
    dNs_sq_dN_par = (N_par*X*Y^2*(1.0 + (Real(mode)*(2.0 - 2*X + (-1.0 + N_par^2)*Y^2))/ (Y^2*sqrtΔ(X, Y, N_par))))/(-1.e0 + X + Y^2)
    return 2.e0 * N .- b .* dNs_sq_dN_par
end


