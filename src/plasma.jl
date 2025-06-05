
struct Plasma{R<:Vector,I<:Extrapolation, T<:Real}
    R_coords::R
    Z_coords::R
    psi_norm_spline::I
    ne_spline::I
    Br_spline::I
    Bz_spline::I
    Bϕ_spline::I
    psi_prof_max::T
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
    psi_norm_spline = CubicSplineInterpolation((r_range, z_range), psi_norm_data; extrapolation_bc=Line())
    psi_range = range(psi_ne[1], psi_ne[end], length(psi_ne))
    ne2 = IMAS.interp1d(psi_ne, ne_prof, :cubic).(psi_range)
    ne_prof_spline = CubicSplineInterpolation(psi_range, log.(ne2); extrapolation_bc=Line())
    ne_data = reshape(ne_prof_spline(reshape(psi_norm_data, length(psi_norm_data))), size(psi_norm_data))
    ne_spline = CubicSplineInterpolation((r_range, z_range), ne_data; extrapolation_bc=Line())
    Br_spline = CubicSplineInterpolation((r_range, z_range), Br_data; extrapolation_bc=Line())
    Bz_spline = CubicSplineInterpolation((r_range, z_range), Bz_data; extrapolation_bc=Line())
    Bϕ_spline = CubicSplineInterpolation((r_range, z_range), Bϕ_data; extrapolation_bc=Line())


    return Plasma(
        R_coords,
        Z_coords,
        psi_norm_spline,
        ne_spline,
        Br_spline,
        Bz_spline,
        Bϕ_spline,
        maximum(psi_ne))
end


function evaluate(spl::Extrapolation, x::AbstractVector{<:Real})
    r = hypot(x[1], x[2])
    z = x[3]
    return spl(r, z)
end


"""
    B_spline(plasma::Plasma, r, z)

Returns total B from the spline representation of individual components
"""
function B_spline(plasma::Plasma, x::AbstractVector{<:Real})
    Br = evaluate(plasma.Br_spline, x)
    Bϕ = evaluate(plasma.Bϕ_spline, x)
    Bz = evaluate(plasma.Bz_spline, x)
    phi = atan(x[2], x[1])
    Bx = Br * cos(phi) - Bϕ * sin(phi)
    By = Br * sin(phi) + Bϕ * cos(phi)
    return vec([Bx By Bz])
end

function n_e(plasma::Plasma, x::AbstractVector{<:Real})
    return exp(evaluate(plasma.ne_spline, x)) #
end
