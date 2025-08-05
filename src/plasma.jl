
struct Plasma{R<:Vector, I2<:Extrapolation, I1<:Extrapolation, T<:Real}
    R_coords::R
    Z_coords::R
    psi_norm_spline::I2
    ne_spline::I2
    Te_spline::I2
    Br_spline::I2
    Bz_spline::I2
    Bϕ_spline::I2
    dvolume_dpsi_spline::I1
    psi_prof_max::T
    
end

function make_2d_prof_spline(r_range, z_range, psi, prof, psi_norm_data)
    psi_range = range(psi[1], psi[end], length(psi))
    prof2 = IMAS.interp1d(psi, prof, :cubic).(psi_range)
    prof_spline = cubic_spline_interpolation(psi_range, log.(prof2); extrapolation_bc=Line())
    prof_data = reshape(prof_spline(reshape(psi_norm_data, length(psi_norm_data))), size(psi_norm_data))
    return cubic_spline_interpolation((r_range, z_range), prof_data; extrapolation_bc=Line())
end

"""
    Plasma(R_coords::Vector{T}, Z_coords::Vector{T}, psi_norm_data::Matrix{T}, psi_ne::Vector{T},
           ne_prof:: Vector{T}, Br_data::Matrix{T}, Bz_data::Matrix{T}, Bϕ_data::Matrix{T}) where {T<:Real}

Generate Plasma structure from 2D maps of densities and magnetic fields
"""
function Plasma(R_coords::Vector{T}, Z_coords::Vector{T}, psi_norm_data::Matrix{T}, psi_prof::Vector{T},
                ne_prof:: Vector{T}, Te_prof:: Vector{T}, Br_data::Matrix{T}, Bz_data::Matrix{T}, 
                Bϕ_data::Matrix{T}, eqt1d_psi:: Vector{T}, eqt1d_psidvolume_dpsi:: Vector{T}) where {T<:Real}
    # Interpolation objects
    r_range = range(R_coords[1], R_coords[end], length(R_coords))
    z_range = range(Z_coords[1], Z_coords[end], length(Z_coords))
    psi_norm_spline = cubic_spline_interpolation((r_range, z_range), psi_norm_data; extrapolation_bc=Line())
    ne_spline = make_2d_prof_spline(r_range, z_range, psi_prof, ne_prof, psi_norm_data)
    Te_spline = make_2d_prof_spline(r_range, z_range, psi_prof, Te_prof, psi_norm_data)
    Br_spline = cubic_spline_interpolation((r_range, z_range), Br_data; extrapolation_bc=Line())
    Bz_spline = cubic_spline_interpolation((r_range, z_range), Bz_data; extrapolation_bc=Line())
    Bϕ_spline = cubic_spline_interpolation((r_range, z_range), Bϕ_data; extrapolation_bc=Line())
    psi_range = range(eqt1d_psi[1], eqt1d_psi[end], length(eqt1d_psi))
    dvolume_dpsi_eq_dist = IMAS.interp1d(eqt1d_psi, eqt1d_psidvolume_dpsi, :cubic).(psi_range)
    dvolume_dpsi_spline = cubic_spline_interpolation(psi_range, dvolume_dpsi_eq_dist; extrapolation_bc=Line())


    return Plasma(
        R_coords,
        Z_coords,
        psi_norm_spline,
        ne_spline,
        Te_spline,
        Br_spline,
        Bz_spline,
        Bϕ_spline,
        dvolume_dpsi_spline,
        maximum(psi_prof))
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

function T_e(plasma::Plasma, x::AbstractVector{<:Real})
    return exp(evaluate(plasma.Te_spline, x)) #
end

function power_deposition_profile(plasma::Plasma, s::Vector{T}, x::Vector{Vector{T}}, dP_ds::Vector{T}, 
                                  psi_dP_dV::Vector{T}, dV_dpsi::Vector{T}) where {T<:Real}
    
    n_points = length(s)
    R = [hypot(x[i][1], x[i][2]) for i in 1:n_points]
    z = [x[i][3] for i in 1:n_points]
    
    psi_s = [plasma.psi_norm_spline(R[i], z[i]) for i in 1:n_points]
    
    # Use Dierckx splines which handle non-uniform spacing better
    dP_ds_spline = Spline1D(s, dP_ds, k=3)  # Cubic spline
    
    dP_dV = zeros(T, length(psi_dP_dV))
    
    for (j, psi_target) in enumerate(psi_dP_dV)
        
        # Create a function that represents psi_spline(s) - psi_target = 0
        # Use Dierckx's root finding capability
        roots_s = T[]
        
        psi_diff_spline = Spline1D(s, psi_s .- psi_target, k=3)
        roots_s = roots(psi_diff_spline)
        power_sum = zero(T)
        for s_root in roots_s
            power_sum += dP_ds_spline(s_root)
        end
        
        dP_dV[j] = power_sum / dV_dpsi[j]
    end
    
    return dP_dV
end


