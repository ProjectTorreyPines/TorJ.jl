"""
    launch_peripheral_rays(x0, N0, width, N_rings, N_azimuthal, dist_focus, focus_shift, one_sigma_width, freq)

Launch peripheral rays for geometrical optics raytracing in a circular pattern.
This function spans a beam perpendicular to the central ray using concentric circles
to discretize the beam cross-section.

# Arguments
- `x0::AbstractVector{T}`: Central ray launch position [x, y, z] in meters
- `N0::AbstractVector{T}`: Central ray direction vector (normalized)
- `width::T`: FWHM beam width in meters
- `N_rings::Integer`: Number of concentric rings
- `N_azimuthal_outer::Integer`: Number of azimuthal points on the outermost ring
- `dist_focus::T`: Distance to focus point in meters
- `focus_shift::T`: Additional focus shift in meters
- `one_sigma_width::Bool`: Whether to use one-sigma width Gaussian weights
- `freq::T`: Frequency in Hz

# Returns
- `ray_positions::Vector{Vector{T}}`: Launch positions for all peripheral rays
- `ray_directions::Vector{Vector{T}}`: Direction vectors for all peripheral rays  
- `ray_weights::Vector{T}`: Statistical weights for each ray
"""
function launch_peripheral_rays(x0::AbstractVector{T}, N0::AbstractVector{T}, 
                                w::T, N_rings::Integer, N_azimuthal_outer::Integer,
                                inverse_curvature_radius::T,  f::T) where {T<:Real}
    
    # Normalized N0 (should already be normalized but better to make sure)
    n0 = N0 ./ LinearAlgebra.norm(N0)
    if isfinite(inverse_curvature_radius)
        R_curv = 1.0 / inverse_curvature_radius 
        λ = constants.c/f
        # beam waist in vacuum - using https://physics.stackexchange.com/questions/270249/for-a-gaussian-beam-how-can-i-calculate-w-0-z-and-z-r-from-knowledge-of
        w0 = (λ * abs(R_curv) * w) / sqrt(λ^2 * R_curv^2 + π^2 * w^4)
        # Distance to f1ocus (not the same as curvature radius!). See stackexchange post above
        ζ_waist = π^2 * R_curv * w^4 /  (λ^2 * R_curv^2 + π^2 * w^4)
        # This respects the sign of the curvature radius
        x_waist = x0 .+ n0 * ζ_waist
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
    r_pts, r_weights = FastGaussQuadrature.gausshermite(N_rings)
    # Integral from +-∞ with e^(-x^2) -> x = √2 r/w
    r_pts .*= sqrt(2.0) * r_pts/w
    # dr/dx = w / √2 | ∫_0^∞ -> 1/2 ∫_{-∞}^∞   
    r_weights .*= w / sqrt(8.0)

    ρ_ring = 2.0*π * max(r_pts)*w / float64(N_azimuthal_outer)
    N_total_rays = 1
    # Need a ragged array to store the points and weights
    for i in 1:N_rings - 1
        N_θ[i] = max(1, round(2π * r_pts[i] * ρ_ring))
        N_total_rays += N_θ[i]
    end
    ray_positions = Vector{Float64}(N_total_rays, 3)
    ray_directions = Vector{Float64}(N_total_rays, 3)
    ray_weights = Vector{Float64}(N_total_rays)
    # Ray index
    k = 0
    for i in 1:N_rings - 1
        θ_azimuthal, θ_weights = FastGaussQuadrature.gausslegendre(N_θ[i])
        θ_azimuthal .*= π
        θ_weights .*= π
        for j in 1:N_θ[i]
            # Convert ray position to Tokamak x-y
            # First calculate local χ and υ in the beam reference frame
            χ = r_weights[i] * cos(θ_azimuthal[j])
            υ = r_weights[i] * sin(θ_azimuthal[j])
            # Calculate ray position in lab frame at launch position
            ray_positions[k + j] .= χ .* e_χ + υ .* e_υ .+ x0
            # Now calculate ray position at the beam waist
            # Use ray direction as temporary array
            if isfinite(inverse_curvature_radius)
                ray_directions[k + j] .= w0/w .* (χ  .* e_χ + υ .* e_υ) .+ x_waist
                if inverse_curvature_radius > 0
                    # Convergent beam -> waist in front of launch position w.r.t. n0
                    ray_directions[k + j] .-= x0
                else
                    # Divergent beam  -> waist behind of launch position w.r.t. n0
                    ray_directions[k + j] .= x0 .- ray_directions[k + j]
                end
            else
                # Perfectly paraxial beam
                ray_directions[k + j] .= n0
            end
            # Calculate complete ray weight from both integrations
            # r comes from the transformation to polar coordinates for this integral
            ray_weights[k + j] = r_pts[i] * r_weights*[i] * θ_weights[j]
        end
        k += N_θ[i]
    end

    return ray_positions, ray_directions, ray_weights
end