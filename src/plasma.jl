
struct Plasma{R<:Vector,I<:Extrapolation}
    R_coords::R
    Z_coords::R
    psi_norm_spline::I
    log_ne_spline::I
    Br_spline::I
    Bz_spline::I
    Bϕ_spline::I
end

function f_interp(itp::Extrapolation, xq::Real, x_max::Real)
    if xq <= x_max
        return itp(xq)
    else
        return 
    end
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
    # Need to smoothly extend n_e
    psi_range = range(psi_ne[1], psi_ne[end], length(psi_ne))
    ne_prof_spline = CubicSplineInterpolation(psi_range, ne_prof; extrapolation_bc=Flat())
    ne_data = reshape(ne_prof_spline(reshape(psi_norm_data, length(psi_norm_data))), size(psi_norm_data))
    dne_dpsi_norm = gradient(ne_prof_spline, psi_ne[end])[1]
    if dne_dpsi_norm > 0.0
        throw(ErrorException("N_e must be decreasing or flat at outermost profile point"))
    end
    ne_lim = 1.e15 # SOL limit
    ne_max = ne_prof[end]
    # Exponential decay
    alpha = -dne_dpsi_norm / (ne_max - ne_lim)
    ne_data[psi_norm_data .>  psi_ne[end]] .= ne_lim .+ (ne_max - ne_lim) .* exp.(-alpha .* (psi_norm_data[psi_norm_data .>  psi_ne[end]] .- psi_ne[end]))
    log_ne_spline = CubicSplineInterpolation((r_range, z_range), log.(ne_data); extrapolation_bc=Flat())
    Br_spline = CubicSplineInterpolation((r_range, z_range), Br_data; extrapolation_bc=Flat())
    Bz_spline = CubicSplineInterpolation((r_range, z_range), Bz_data; extrapolation_bc=Flat())
    Bϕ_spline = CubicSplineInterpolation((r_range, z_range), Bϕ_data; extrapolation_bc=Flat())

    
    return Plasma(
        R_coords,
        Z_coords,
        psi_norm_spline,
        log_ne_spline,
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
