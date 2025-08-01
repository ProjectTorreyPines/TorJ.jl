module TorJPlotsExt

using TorJ
using Plots
using PyPlot
using TorJ: LinearAlgebra
using TorJ: IMAS
using Interpolations

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
    steering_angle_pol = deg2rad(30.0) # Convert from TORBEAM convetion to IMAS (they are the same for phi_tor ==0)
    steering_angle_tor = 0.0  # Radial injection
    N0 = collect(IMAS.pol_tor_angles_2_vector(steering_angle_pol, steering_angle_tor))  # Direction vector
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
    
    # Use PyPlot backend for better window display
    pyplot()
    
    # First try a static plot to test if window opens
    p = plot3d(xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
               title="Peripheral Rays ($N_rings rings)",
               legend=:topright, dpi=300, aspect_ratio=:equal; kwargs...)
    
    # Normalize weights for alpha (0.1 to 1.0 range)
    min_weight, max_weight = extrema(ray_weights)
    normalized_weights = (ray_weights .- min_weight) ./ (max_weight - min_weight)
    alpha_values = 0.05 .+ 0.95 .* normalized_weights  # Scale to 0.1-1.0 range
    
    # Plot each ray as a line segment with color and alpha based on weight
    for i in 1:N_rays
        weight_color = normalized_weights[i]
        alpha_val = alpha_values[i]
        plot3d!(p, [x_origins[i], x_ends[i]], [y_origins[i], y_ends[i]], [z_origins[i], z_ends[i]],
                color=:viridis, line_z=weight_color, linewidth=2, alpha=alpha_val, label="")
    end
    
    # Add origin points as scatter
    scatter3d!(p, x_origins, y_origins, z_origins,
            zcolor=normalized_weights, c=:viridis, markersize=4, markerstrokewidth=0,
            label="Ray origins", colorbar_title="Ray Weight", camera=(45, 30))
    
    # Display the static plot
    display(p)
    
    return p 
    
    # Print statistics
    println("3D Ray Plot Statistics:")
    println("Total rays: $(N_rays)")
    println("Ray weights - Min: $(minimum(ray_weights)), Max: $(maximum(ray_weights))")
    println("Sum of all weights: $(sum(ray_weights))")
    
end

"""
    plot_beam_trajectories_3d(trajectories, weights; kwargs...)

Plot beam trajectories in 3D with color coding based on ray weights.

# Arguments
- `trajectories`: Vector of ray trajectory data from make_beam
- `weights`: Ray weights corresponding to each trajectory
- `kwargs...`: Additional arguments passed to plot

# Returns
- Plots.jl 3D plot object
"""
function plot_beam_trajectories_3d(arc_lengths, trajectories, ray_powers, weights; kwargs...)
    # Use PyPlot backend for better window display
    pyplot()
    
    # Create the 3D plot
    p = plot3d(xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
               title="Beam Trajectories with Weight Color Coding",
               legend=false, dpi=300, aspect_ratio=:equal; kwargs...)
    println("Got $(length(weights)) rays")
    total_points_before = sum(length(ray_trajectory) for ray_trajectory in trajectories)
    println("Total trajectory points before resampling: $(total_points_before)")
    
    # Resample trajectories to 0.01 m resolution
    ds_target = 0.01  # Target spacing in meters
    resampled_trajectories = []
    total_points_after = 0
    
    for i in 1:length(weights)
        arc_length = arc_lengths[i]
        ray_trajectory = trajectories[i]
        
        println("Ray $(i): $(length(ray_trajectory)) points before resampling")
        
        if length(ray_trajectory) < 2 || length(arc_length) < 2
            # Handle edge case of very short trajectories
            push!(resampled_trajectories, ray_trajectory)
            total_points_after += length(ray_trajectory)
            println("Ray $(i): $(length(ray_trajectory)) points after resampling (too short to resample)")
            continue
        end
        
        # Create interpolation functions for x, y, z coordinates
        ray_x = [point[1] for point in ray_trajectory]
        ray_y = [point[2] for point in ray_trajectory]  
        ray_z = [point[3] for point in ray_trajectory]
        
        # Create new arclength grid with target spacing
        s_min = arc_length[1]
        s_max = arc_length[end]
        n_new_points = max(2, Int(ceil((s_max - s_min) / ds_target)) + 1)
        s_new = range(s_min, s_max, length=n_new_points)
        
        # Interpolate coordinates
        println(arc_length)
        interp_x = Interpolations.linear_interpolation(arc_length, ray_x)
        interp_y = Interpolations.linear_interpolation(arc_length, ray_y)
        interp_z = Interpolations.linear_interpolation(arc_length, ray_z)
        
        # Generate resampled trajectory
        resampled_trajectory = [[interp_x(s), interp_y(s), interp_z(s)] for s in s_new]
        push!(resampled_trajectories, resampled_trajectory)
        
        total_points_after += length(resampled_trajectory)
        println("Ray $(i): $(length(resampled_trajectory)) points after resampling")
    end
    
    println("Total trajectory points after resampling: $(total_points_after)")
    
    # Normalize weights for color mapping
    min_weight, max_weight = extrema(weights)
    normalized_weights = (weights .- min_weight) ./ (max_weight - min_weight)
    alpha_values = 0.05 .+ 0.95 .* normalized_weights  # Scale to 0.1-1.0 range
    
    # Plot each resampled trajectory
    for i in 1:length(weights)
        println("Initial vs. final beam power: $(ray_powers[i][1]) $(ray_powers[i][end])")
        
        # Use resampled trajectory
        ray_trajectory = resampled_trajectories[i]
        
        # Extract coordinates from resampled trajectory points
        ray_x = [point[1] for point in ray_trajectory]
        ray_y = [point[2] for point in ray_trajectory]  
        ray_z = [point[3] for point in ray_trajectory]
        
        # Color based on normalized weight
        plot3d!(p, ray_x, ray_y, ray_z,
               color=:viridis, line_z=fill(normalized_weights[i], length(ray_x)),
               linewidth=2, alpha=alpha_values[i], label="")
    end
    
    # Add colorbar
    scatter3d!(p, [0], [0], [0], zcolor=[0], c=:viridis, 
              markersize=0, markerstrokewidth=0, colorbar_title="Ray Weight")
    
    # Display the plot
    display(p)
    
    return p
end

"""
    plot_beam_from_setup(; s_max=0.4, kwargs...)

Create and plot beam trajectories using parameters from setup.jl.

# Arguments  
- `s_max`: Maximum integration distance (default: 0.4)
- `kwargs...`: Additional arguments passed to plot

# Returns
- Plots.jl 3D plot object
"""
function plot_beam_from_setup(; s_max=0.4, kwargs...)
    # Import setup parameters (similar to test_make_beam.jl)
    include(joinpath(@__DIR__, "../test/tests/setup.jl"))
    
    # Generate beam trajectories using setup parameters
    arc_lengths, trajectories, ray_powers, ray_weights = TorJ.make_beam(plasma, R0, phi0, z0, steering_angle_tor,
                                               steering_angle_pol, spot_size, 
                                               inverse_curvature_radius, freq, 1, s_max)
    
    # Call the main plotting function
    return plot_beam_trajectories_3d(arc_lengths, trajectories, ray_powers, ray_weights; kwargs...)
end

end # module TorJPlotsExt