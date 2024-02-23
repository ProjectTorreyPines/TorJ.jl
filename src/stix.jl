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

"""
    ω_pe(ne::Real)

Returns electron plasma frequency [rad/s] given electron density in m⁻³
"""
function ω_pe(ne::Real)
    return sqrt(ne * constants.e^2 / (constants.ϵ_0 * constants.m_e))
end

"""
    ω_ce(B::Real)

Returns electron cyclotron frequency [rad/s] given magnetic field B in T
"""
function ω_ce(B::Real)
    return constants.e * abs(B) / constants.m_e
end

"""
    B_ω_ce(ω::Real)

Returns magnetic field B in T for a given electron cyclotron frequency [rad/s]
"""
function B_ω_ce(ω::Real)
    return ω / constants.e * constants.m_e
end

"""
    ω_ci(B::Real, Z::Real, A::Real)

Returns ion cyclotron frequency [rad/s] given magnetic field B in T and the ion charge and mass in amu
"""
function ω_ci(B::Real, Z::Real, A::Real)
    return constants.e * abs(B) * Z / A * constants.m_p
end

"""
    S_stix(plasma::Plasma, r, z, ω)

S (Sum term) of the Stix dielectric tensor
"""
function S_stix(plasma::Plasma, r, z, ω)
    ne = plasma.ne_spline(r, z)
    B = B_spline(plasma, r, z)
    ωpe = ω_pe(ne)
    ωce = ω_ce(B)
    return 1.0 - ωpe^2 / (ω^2 - ωce^2)
end

"""
    D_stix(plasma::Plasma, r, z, ω)

D (Difference term) of the Stix dielectric tensor
"""
function D_stix(plasma::Plasma, r, z, ω)
    ne = plasma.ne_spline(r, z)
    B = B_spline(plasma, r, z)
    ωpe = ω_pe(ne)
    ωce = ω_ce(B)
    return ωpe^2 * ωce / (ω * (ω^2 - ωce^2))
end

"""
    P_stix(plasma::Plasma, r, z, ω)

P (Plasma term) of the Stix dielectric tensor
"""
function P_stix(plasma::Plasma, r, z, ω)
    ne = plasma.ne_spline(r, z)
    ωpe = ω_pe(ne)
    return 1.0 - ωpe^2 / ω^2
end