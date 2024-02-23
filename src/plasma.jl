struct Plasma{R<:AbstractRange,I<:AbstractInterpolation}
    R_coords::R
    Z_coords::R
    ne_spline::I
    Br_spline::I
    Bz_spline::I
    Bϕ_spline::I
end

"""
    Plasma(R_coords::AbstractRange, Z_coords::AbstractRange, ne_data::Matrix{T}, Br_data::Matrix{T}, Bz_data::Matrix{T}, Bϕ_data::Matrix{T})

Generate Plasma structure from 2D maps of densities and magnetic fields
"""
function Plasma(R_coords::AbstractRange, Z_coords::AbstractRange, ne_data::Matrix{T}, Br_data::Matrix{T}, Bz_data::Matrix{T}, Bϕ_data::Matrix{T}) where {T<:Real}
    # Interpolation objects
    ne_spline = CubicSplineInterpolation((R_coords, Z_coords), ne_data; extrapolation_bc=Flat())
    Br_spline = CubicSplineInterpolation((R_coords, Z_coords), Br_data; extrapolation_bc=Flat())
    Bz_spline = CubicSplineInterpolation((R_coords, Z_coords), Bz_data; extrapolation_bc=Flat())
    Bϕ_spline = CubicSplineInterpolation((R_coords, Z_coords), Bϕ_data; extrapolation_bc=Flat())

    return Plasma(
        R_coords,
        Z_coords,
        ne_spline,
        Br_spline,
        Bz_spline,
        Bϕ_spline)
end

"""
    B_spline(plasma::Plasma, r, z)

Returns total B from the spline representation of individual components
"""
function B_spline(plasma::Plasma, r, z)
    Br = plasma.Br_spline(r, z)
    Bϕ = plasma.Bϕ_spline(r, z)
    Bz = plasma.Bz_spline(r, z)
    B = sqrt(Br^2 + Bϕ^2 + Bz^2)
    return B
end