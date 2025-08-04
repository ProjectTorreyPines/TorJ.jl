"""
    launch_peripheral_rays(x0, N0, width, N_rings, N_azimuthal, dist_focus, focus_shift, one_sigma_width, f)

Launch peripheral rays for geometrical optics raytracing in a circular pattern.
This function spans a beam perpendicular to the central ray using concentric circles
to discretize the beam cross-section.

# Arguments (All use IMAS conventions where applicable)
- `x0::AbstractVector{T}`: Central ray launch position [x, y, z] in meters
- `N0::AbstractVector{T}`: Central ray direction vector (normalized)
- `w::T`: beam width in meters
- `inverse_curvature_radius`: Inverse curvature radius of the circular beam at launch
- `N_rings::Integer`: Number of concentric rings
- `f::T`: Frequency in Hz
- `min_azimuthal_points::Integer`: Minimal number of azimuthal points should be at least 5

# Returns
- `ray_positions::Vector{Vector{T}}`: Launch positions for all peripheral rays
- `ray_directions::Vector{Vector{T}}`: Direction vectors for all peripheral rays  
- `ray_weights::Vector{T}`: Statistical weights for each ray
- `normalize_weight_sum`: Make sure all weights sum up to one. Should only be false for tests.
"""
function launch_peripheral_rays(x0::AbstractVector{T}, N0::AbstractVector{T}, 
                                w::T, N_rings::Integer, inverse_curvature_radius::T,  
                                f::T; min_azimuthal_points = 5, normalize_weight_sum=true) where {T<:Real}
     if N_rings < 2
        throw(ArgumentError("N_rings = $N_rings < 2 which is the minimum"))
    end
    # Normalized N0 (should already be normalized but better to make sure)
    n0 = N0 ./ LinearAlgebra.norm(N0)
    # Define a minimum of azimuthal points
    
    if isfinite(inverse_curvature_radius)
        R_curv = 1.0 / inverse_curvature_radius 
        λ = constants.c/f
        # beam waist in vacuum - using https://physics.stackexchange.com/questions/270249/for-a-gaussian-beam-how-can-i-calculate-w-0-z-and-z-r-from-knowledge-of
        w0 = (λ * abs(R_curv) * w) / sqrt(λ^2 * R_curv^2 + π^2 * w^4)
        # Distance to waist (not the same as curvature radius!). See stackexchange post above
        z_waist = π^2 * R_curv * w^4 /  (λ^2 * R_curv^2 + π^2 * w^4)
        # This respects the sign of the curvature radius
        x_waist = x0 .+ n0 .* z_waist
    else
        # Completely paraxial beam
        w0 = w
    end
    
    # First vector perpendicular to N0
    e_χ = zeros(T, 3)
    e_χ[1] = 1.0
    e_χ[2] = 0.0
    e_χ[3] = -n0[1] / n0[3]
    
    # Second vector perpendicular to N0
    # The χ-υ basis is used for transformation to Cartesian lab frame
    e_υ = zeros(T, 3)
    e_υ[1] = -n0[1] * n0[2] / n0[3]
    e_υ[2] = n0[3] - n0[1]
    e_υ[3] = -n0[2]
    
    # Normalize orthogonal vectors
    e_χ ./= LinearAlgebra.norm(e_χ)
    e_υ ./= LinearAlgebra.norm(e_υ)
    
    # Create points for the integration
    # Remove the first positive r since it is close to zero
    r_pts, r_weights = (v[N_rings+2:end] for v in FastGaussQuadrature.gausshermite(N_rings*2+2)) # use even integer
    # Integral from +-∞ with e^(-x^2) -> r = w/√2 * x
    r_pts .*= w/sqrt(2.0)
    # dr/dx = w / √2
    r_weights .*= w / sqrt(2.0)
    N_θ = zeros(Int64, N_rings)
    N_total_rays = Int64(0)
    # Need a ragged array to store the points and weights
    for i in 1:N_rings
        N_θ[i] = max(1, round(Int64, min_azimuthal_points * r_pts[i] / r_pts[1]))
        N_total_rays += N_θ[i]
    end
    ray_positions = zeros(Float64, N_total_rays, 3)
    ray_directions = zeros(Float64, N_total_rays, 3)
    ray_weights = zeros(Float64, N_total_rays)
    # Ray index
    k = 0
    for i in 1:N_rings
        # Simple trapezoid here since we are almost constant in θ
        θ_azimuthal = [2π * Float64(l) / Float64(N_θ[i]) for l in 0:(N_θ[i]-1)]
        θ_weights = fill(2π /Float64(N_θ[i]), N_θ[i])
        for j in 1:N_θ[i]
            # Convert ray position to Tokamak x-y
            # First calculate local χ and υ in the beam reference frame
            χ = r_pts[i] * cos(θ_azimuthal[j])
            υ = r_pts[i] * sin(θ_azimuthal[j])
            # Calculate ray position in lab frame at launch position
            ray_positions[k + j,:] .= χ .* e_χ .+ υ .* e_υ .+ x0
            # Now calculate ray position at the beam waist
            # Use ray direction as temporary array
            if isfinite(inverse_curvature_radius)
                ray_directions[k + j,:] .= w0/w .* (χ  .* e_χ + υ .* e_υ) .* sign(inverse_curvature_radius)
                if inverse_curvature_radius > 0
                    # Convergent beam -> waist in front of launch position w.r.t. n0
                    ray_directions[k + j,:] .+= x_waist - ray_positions[k + j,:]
                else
                    # Divergent beam  -> waist behind of launch position w.r.t. n0
                    ray_directions[k + j,:] .+= ray_positions[k + j,:] .- x_waist
                end
                ray_directions[k + j,:] ./= LinearAlgebra.norm(ray_directions[k + j,:])
            else
                # Perfectly paraxial beam
                ray_directions[k + j,:] .= n0
            end
            # Calculate complete ray weight from both integrations
            # r comes from the transformation to polar coordinates for this integral
            ray_weights[k + j] = r_pts[i] * r_weights[i] * θ_weights[j]
        end
        k += N_θ[i]
    end
    # Add normalization
    if normalize_weight_sum
        ray_weights ./= sum(ray_weights)
    else
        ray_weights .*= 2.0 / (w^2* π)
    end
    
    return ray_positions, ray_directions, ray_weights
end