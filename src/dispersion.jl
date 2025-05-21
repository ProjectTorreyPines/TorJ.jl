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
    B = B_spline(plasma, r, phi, x[3])
    B_abs = LinearAlgebra.norm(B)
    n_e = plasma.ne_spline(r, x[3])
    N_par = LinearAlgebra.dot(N, B ./ B_abs) 
    X = n_e * constants.e^2.e0/(constants.ϵ_0 * constants.m_e * omega^2)
    Y = constants.e / (constants.m_e * B_abs * omega)
    return X, Y, N_par
end

function refractive_index_sq( X::Real, Y::Real, N_par::Real, mode::Integer)
    Δ = (1.0 - N_par^2)^2 + 4.0 * N_par^2 * (1.0 - X) / Y^2
    if Δ < 0
        return 0
    end
    Δ = sqrt(Δ)
    Ns_sq = 1.e0 - X + (1.0 + Real(mode) * Δ + N_par^2)/(2.0 * (-1.0 + X + Y^2)) * X * Y^2
    return Ns_sq
end

function dispersion_relation(x::AbstractVector{<:Real}, N::AbstractVector{<:Real}, plasma::Plasma, omega:: Real, mode::Integer)
    N_abs = LinearAlgebra.norm(N)
    X, Y, N_par = eval_plasma(plasma, x, N, omega)
    Ns_sq = refractive_index_sq(X, Y, N_par, mode)
    return N_abs^2 - Ns_sq
end

