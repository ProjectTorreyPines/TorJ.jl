
struct Plasma{R<:Vector, I2<:Extrapolation, I1<:Extrapolation, T<:Real}
    R_coords::R
    Z_coords::R
    psi_norm_spline::I2
    ne_spline::I2
    Te_spline::I2
    Br_spline::I2
    Bz_spline::I2
    Bϕ_spline::I2
    volume_psi_spline::I1
    psi_prof_max::T
    
end

function make_2d_prof_spline(r_range, z_range, psi_norm, prof, psi_norm_grid)
    psi_norm_range = range(psi_norm[1], psi_norm[end], length(psi_norm))
    prof2 = IMAS.interp1d(psi_norm, prof, :cubic).(psi_norm_range)
    prof_spline = cubic_spline_interpolation(psi_norm_range, log.(prof2); extrapolation_bc=Line())
    prof_data = reshape(prof_spline(reshape(psi_norm_grid, length(psi_norm_grid))), size(psi_norm_grid))
    return cubic_spline_interpolation((r_range, z_range), prof_data; extrapolation_bc=Line())
end

"""
    Plasma(dd::IMAS.dd)

Generate Plasma structure from dd.
"""
function Plasma(dd::IMAS.dd; ne_scale::Float64=1.0)
    eqt = dd.equilibrium.time_slice[]
    eqglobs = eqt.global_quantities
    eqt1d = eqt.profiles_1d
    eqt2d = IMAS.findfirst(:rectangular, eqt.profiles_2d)
    cp1d = dd.core_profiles.profiles_1d[]
    return Plasma(eqt2d.grid.dim1, eqt2d.grid.dim2, eqt2d.psi, 
                  eqt2d.b_field_r, eqt2d.b_field_z, eqt2d.b_field_tor,
                  eqglobs.psi_axis, eqglobs.psi_boundary,
                  eqt1d.psi, eqt1d.volume, 
                  cp1d.grid.psi, cp1d.electrons.density .* ne_scale,
                  cp1d.electrons.temperature)
end

"""
    Plasma(R_coords::Vector{T}, Z_coords::Vector{T}, psi_norm_data::Matrix{T}, psi_ne::Vector{T},
           ne_prof:: Vector{T}, Br_data::Matrix{T}, Bz_data::Matrix{T}, Bϕ_data::Matrix{T}) where {T<:Real}

Generate Plasma structure from 2D maps of densities and magnetic fields. All fields use the IMAS conventions.
"""
function Plasma(R_coords::Vector{T}, Z_coords::Vector{T}, psi_grid::Matrix{T}, 
                Br_data::Matrix{T}, Bz_data::Matrix{T}, Bϕ_data::Matrix{T},
                psi_axis::T, psi_boundary::T,
                eqt1d_psi:: Vector{T}, eqt1d_volume:: Vector{T},
                psi_core_prof::Vector{T}, ne_prof:: Vector{T}, Te_prof:: Vector{T}) where {T<:Real}
    # Interpolation objects
    _norm = psi -> (psi .- psi_axis) ./ (psi_boundary - psi_axis)
    r_range = range(R_coords[1], R_coords[end], length(R_coords))
    z_range = range(Z_coords[1], Z_coords[end], length(Z_coords))
    psi_norm_spline = cubic_spline_interpolation((r_range, z_range), _norm(psi_grid); extrapolation_bc=Line())

    Br_spline = cubic_spline_interpolation((r_range, z_range), Br_data; extrapolation_bc=Line())
    Bz_spline = cubic_spline_interpolation((r_range, z_range), Bz_data; extrapolation_bc=Line())
    Bϕ_spline = cubic_spline_interpolation((r_range, z_range), Bϕ_data; extrapolation_bc=Line())
    
    psi_norm_range = range(_norm(eqt1d_psi[1]), _norm(eqt1d_psi[end]), length(eqt1d_psi))
    volume_eq_dist = IMAS.interp1d(_norm(eqt1d_psi), eqt1d_volume, :cubic).(psi_norm_range)
    volume_psi_spline = cubic_spline_interpolation(psi_norm_range, volume_eq_dist; extrapolation_bc=Line())

    ne_spline = make_2d_prof_spline(r_range, z_range, _norm(psi_core_prof), ne_prof, _norm(psi_grid))
    Te_spline = make_2d_prof_spline(r_range, z_range, _norm(psi_core_prof), Te_prof, _norm(psi_grid))
    
    return Plasma(
        R_coords,
        Z_coords,
        psi_norm_spline,
        ne_spline,
        Te_spline,
        Br_spline,
        Bz_spline,
        Bϕ_spline,
        volume_psi_spline,
        maximum(_norm(eqt1d_psi)))
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
                                  psi_dP_dV::Vector{T}; verbose=false) where {T<:Real}
    
    n_points = length(s)
    R = [hypot(x[i][1], x[i][2]) for i in 1:n_points]
    z = [x[i][3] for i in 1:n_points]
    
    psi_s = [plasma.psi_norm_spline(R[i], z[i]) for i in 1:n_points]
    
    # Use Dierckx splines which handle non-uniform spacing better
    dP_ds_spline = Dierckx.Spline1D(s, dP_ds, k=3)  # Cubic spline
    
    dP_dV = zeros(T, length(psi_dP_dV))
    P = 0.0
    j =  length(psi_dP_dV)
    psi_diff_spline = Dierckx.Spline1D(s, psi_s .- psi_dP_dV[j], k=3)
    # Find root with outer flux surface shell
    outer_roots = Dierckx.roots(psi_diff_spline)
    outer_volume = plasma.volume_psi_spline(psi_dP_dV[j])
    j -= 1
    while j > 0
        # Step through all volumes and integrate the power and volume change for each psi
        inner_volume = plasma.volume_psi_spline(psi_dP_dV[j])
        δV = outer_volume - inner_volume
        # Create a function that represents psi_spline(s) - psi_target = 0
        psi_diff_spline = Dierckx.Spline1D(s, psi_s .- psi_dP_dV[j], k=3)
        # First find intersection with inner flux surface shell
        inner_roots = Dierckx.roots(psi_diff_spline)
        intervals = sort(append!(outer_roots, inner_roots))
        if length(intervals) < 2
            # length(intervals) = 0: No intersection () -> The beam did not go into the flux surface
            # length(intervals) = 1: End of the ray, does not get counted
            # Since we are moving from the outside to the inside once we do not find a intersection we are done
            break
        elseif length(intervals) % 2 != 0
            # length(intervals) % 2 = 1: End of the ray, does not count
            intervals = intervals[1:end-1]
        end
        k = 1
        δP = 0.0
        while k < length(intervals)
            # Integrate all power deposited in this flux surface shell
            # This should be an average but we are expressing <P> averaged over δs which is multiplied by δs/δV, hence the δs cancels
            δP += abs(Dierckx.integrate(dP_ds_spline, intervals[k], intervals[k+1]))
            k += 2
        end
        if verbose
            println("j: $j | Ψ: $(psi_dP_dV[j]) | Intervals: $intervals | δP $δP")
        end
        # Form the differential
        dP_dV[j] = δP/δV
        # This should be very close to 1 - sum_i P_i(s_max) (1 - the power remaining in all rays) if we did it right
        P += δP

        # Prepare the next step
        j -= 1
        outer_volume = inner_volume
        outer_roots = inner_roots
    end
    return dP_dV, P
end


