"""
    dispersion_relation(plasma::Plasma, r, ϕ, z, kr, nϕ, kz, ω)

Cold plasma dispersion relation
"""
function dispersion_relation(plasma::Plasma, r, ϕ, z, kr, nϕ, kz, ω)
    c = constants.c

    kpar2, kper2 = k_par2_per2(plasma, r, ϕ, z, kr, nϕ, kz)

    npar2 = (kpar2 * c^2) / ω^2
    nper2 = (kper2 * c^2) / ω^2

    S = S_stix(plasma, r, z, ω)
    D = D_stix(plasma, r, z, ω)
    P = P_stix(plasma, z, z, ω)

    C4 = S
    C2 = (npar2 - S) * (P + S) + D^2
    C0 = P * ((npar2 - S)^2 - D^2)

    return C4 * nper2^2 + C2 * nper2 + C0
end

"""
    k_par2_per2(plasma::Plasma, r, ϕ, z, kr, nϕ, kz)

Compute kpar² and kper²
"""
function k_par2_per2(plasma::Plasma, r, ϕ, z, kr, nϕ, kz)
    Br = plasma.Br_spline(r, z)
    Bϕ = plasma.Bϕ_spline(r, z)
    Bz = plasma.Bz_spline(r, z)
    B = sqrt(Br^2 + Bϕ^2 + Bz^2)

    kϕ = nϕ / r
    kpar2 = ((kr * Br + kz * Bz + kϕ * Bϕ) / B)^2.0

    k2 = kr^2 + kz^2 + kϕ^2
    kper2 = k2 - kpar2

    return (kpar2=kpar2, kper2=kper2)
end