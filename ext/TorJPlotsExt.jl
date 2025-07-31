module TorJPlotsExt

using TorJ
using Plots
using PyPlot
using TorJ: LinearAlgebra
using TorJ: IMAS

"""
    plot_peripheral_rays_3d(; N_rings=5, N_azimuthal_outer=8, ray_length=1.0, kwargs...)

Generate 3D plot of peripheral rays with weights indicated by color.

# Arguments
- `N_rings`: Number of concentric rings (default: 5)
- `N_azimuthal_outer`: Number of azimuthal points on outermost ring (default: 25)
- `ray_length`: Length of ray segments to plot in meters (default: 1.0)
- `kwargs...`: Additional arguments passed to plot

# Returns
- Plots.jl 3D plot object
"""
function plot_peripheral_rays_3d(; N_rings=5, ray_length=1.0, kwargs...)
    # Parameters from setup.jl and user specification
    x0 = [2.5, 0.0, 0.4]
    α = deg2rad(30.0) 
    β = 0.0               # Radial injection
    N0 = collect(IMAS.pol_tor_angles_2_vector(α, β))  # Direction vector
    w = 0.0174            # Beam width parameter
    inverse_curvature_radius = 1.0/3.99  # Curvature parameter
    freq = 110.0e9         # Frequency
    
    # Generate peripheral rays using the fixed function
    ray_origins, ray_directions, ray_weights = TorJ.launch_peripheral_rays(
        x0, N0, w, N_rings, inverse_curvature_radius, freq)
    
    # Extract data from matrices (N_rays × 3 format)
    N_rays = size(ray_origins, 1)
    x_origins = ray_origins[:, 1]
    y_origins = ray_origins[:, 2]  
    z_origins = ray_origins[:, 3]
    
    # ray_directions now contains normalized direction vectors
    # Calculate end points (ray_length meters along the direction)
    x_ends = x_origins + ray_length * ray_directions[:, 1]
    y_ends = y_origins + ray_length * ray_directions[:, 2] 
    z_ends = z_origins + ray_length * ray_directions[:, 3]
    gr()
    anim = Animation()
    for angle in range(0, stop = 180, step = 1)
        # Create 3D plot with rotating camera
        p = plot3d(xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
                title="Peripheral Rays ($N_rings rings)",
                legend=:topright, dpi=300; kwargs...)
        
        # Plot each ray as a line segment with color based on weight
        for i in 1:N_rays
            weight_color = ray_weights[i]
            plot3d!(p, [x_origins[i], x_ends[i]], [y_origins[i], y_ends[i]], [z_origins[i], z_ends[i]],
                    color=:viridis, line_z=weight_color, linewidth=2, alpha=0.8, label="")
        end
        
        # Add origin points as scatter
        scatter3d!(p, x_origins, y_origins, z_origins,
                zcolor=ray_weights, c=:viridis, markersize=4, markerstrokewidth=0,
                label="Ray origins", colorbar_title="Ray Weight", camera=(angle, 30))
        frame(anim, p)
    end 
    
    # Display the animation directly
    anim
    
    # Print statistics
    println("3D Ray Plot Statistics:")
    println("Total rays: $(N_rays)")
    println("Ray weights - Min: $(minimum(ray_weights)), Max: $(maximum(ray_weights))")
    println("Sum of all weights: $(sum(ray_weights))")
    
    return anim
end

end # module TorJPlotsExt