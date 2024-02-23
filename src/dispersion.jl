"""
    dispersion_relation(plasma::Plasma, r, ϕ, z, kr, nϕ, kz, ω)

Cold plasma dispersion relation
"""
function dispersion_relation(plasma::Plasma, r, ϕ, z, kr, nϕ, kz, ω)
    c = constants.c

    kpar2, kperp2 = k_par2_perp2(plasma, r, ϕ, z, kr, nϕ, kz)

    S = S_stix(plasma, r, z, ω)
    D = D_stix(plasma, r, z, ω)
    P = P_stix(plasma, z, z, ω)

    term1 = S * (kperp2^2 * c^4) / ω^4
    term2 = ((S + P) * (S - (kpar2 * c^2) / ω^2) - D^2) * (kperp2 * c^2) / ω^2
    term3 = P * ((S - (kpar2 * c^2) / ω^2)^2 - D^2)

    return term1 - term2 + term3
end

"""
    k_par2_perp2(plasma::Plasma, r, ϕ, z, kr, nϕ, kz)

Compute kpar² and kperp² 
"""
function k_par2_perp2(plasma::Plasma, r, ϕ, z, kr, nϕ, kz)
    Br = plasma.Br_spline(r, z)
    Bϕ = plasma.Bϕ_spline(r, z)
    Bz = plasma.Bz_spline(r, z)
    B = sqrt(Br^2 + Bϕ^2 + Bz^2)
    kϕ = nϕ / r
    kpar2 = ((kr * Br + kz * Bz + kϕ * Bϕ) / B)^2.0

    k2 = kr^2 + kz^2 + kϕ^2
    kperp2 = k2 - kpar2

    return (kpar2=kpar2, kperp2=kperp2)
end