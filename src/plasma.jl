
struct Plasma{R<:Vector,I<:Extrapolation}
    R_coords::R
    Z_coords::R
    psi_norm_spline::I
    ne_spline::I
    Br_spline::I
    Bz_spline::I
    Bϕ_spline::I
end

"""
    Plasma(R_coords::Vector{T}, Z_coords::Vector{T}, psi_norm_data::Matrix{T}, psi_ne::Vector{T},
           ne_prof:: Vector{T}, Br_data::Matrix{T}, Bz_data::Matrix{T}, Bϕ_data::Matrix{T}) where {T<:Real}

Generate Plasma structure from 2D maps of densities and magnetic fields
"""
function Plasma(R_coords::Vector{T}, Z_coords::Vector{T}, psi_norm_data::Matrix{T}, psi_ne::Vector{T},
                ne_prof:: Vector{T}, Br_data::Matrix{T}, Bz_data::Matrix{T}, Bϕ_data::Matrix{T}) where {T<:Real}
    # Interpolation objects
    r_range = range(R_coords[1], R_coords[end], length(R_coords))
    z_range = range(Z_coords[1], Z_coords[end], length(Z_coords))
    psi_norm_spline = CubicSplineInterpolation((r_range, z_range), psi_norm_data; extrapolation_bc=Flat())
    psi_range = range(psi_ne[1], psi_ne[end], length(psi_ne))
    ne_prof_spline = CubicSplineInterpolation(psi_range, ne_prof; extrapolation_bc=Flat())
    ne_data = reshape(ne_prof_spline(reshape(psi_norm_data, length(psi_norm_data))), size(psi_norm_data))
    ne_spline = CubicSplineInterpolation((r_range, z_range), ne_data; extrapolation_bc=Flat())
    Br_spline = CubicSplineInterpolation((r_range, z_range), Br_data; extrapolation_bc=Flat())
    Bz_spline = CubicSplineInterpolation((r_range, z_range), Bz_data; extrapolation_bc=Flat())
    Bϕ_spline = CubicSplineInterpolation((r_range, z_range), Bϕ_data; extrapolation_bc=Flat())

    
    return Plasma(
        R_coords,
        Z_coords,
        psi_norm_spline,
        ne_spline,
        Br_spline,
        Bz_spline,
        Bϕ_spline)
end

"""
    B_spline(plasma::Plasma, r, z)

Returns total B from the spline representation of individual components
"""
function B_spline(plasma::Plasma, r, phi, z)
    Br = plasma.Br_spline(r, z)
    Bϕ = plasma.Bϕ_spline(r, z)
    Bz = plasma.Bz_spline(r, z)
    Bx = Br * cos(phi) - Bϕ * sin(phi)
    By = Br * sin(phi) + Bϕ * cos(phi)
    return vec([Bx By Bz])
end
