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

function eval_plasma(plasma::Plasma, x::T, y::T, z::T, Nx::T, Ny::T, Nz::T, omega::T) where {T<:Real}
    N_abs = sqrt(Nx^2 + Ny^2 + Nz^2)
    r = sqrt(x^2 + y^2)
    phi = atan(y,x)
    Bx, By, Bz = B_spline(plasma, r, phi, z)
    n_e = plasma.ne_spline(r, z)
    B_abs = sqrt(Bx^2 + By^2 + Bz^2)
    N_par = Bx*Nx + By*Ny + Bz*Nz
    N_par /= B_abs * N_abs
    X = n_e * constants.e^2.e0/(constants.ϵ_0 * constants.m_e * omega^2)
    Y = constants.e / (constants.m_e * B_abs * omega)
    return X, Y, N_par
end

function refractive_index_sq( X::T, Y::T, N_par::T, mode::Integer) where {T<:Real}
    Δ = (1.0 - N_par^2)^2 + 4.0 * N_par^2 * (1.0 - X) / Y^2
    Δ = sqrt(Complex(Δ))
    Ns_sq = 1.e0 - X + (1.0 + Real(mode) * Δ + N_par^2)/(2.0 * (-1.0 + X + Y^2)) * X * Y^2
    return Ns_sq
end

function dispersion_relation(x::T, N::T, plasma::Plasma, omega:: T, mode::Integer) where {T<:Real}
    N_abs = sqrt(N[1]^2 + N[2]^2 + N[2]^2)
    X, Y, N_par = eval_plasma(plasma, x[1], x[2], x[3], N[1], N[2], N[3], omega)
    Ns_sq = refractive_index_sq(X, Y, N_par, mode)
    return N_abs^2 - Ns_sq
end

