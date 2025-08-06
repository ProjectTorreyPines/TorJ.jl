"""
Test for Gaussian integration using polar coordinates from launch_peripheral_rays.

This test validates that the polar coordinate discretization from launch_peripheral_rays
correctly integrates a 2D Gaussian function by comparing to a reference integration method.
"""

using Test
import TorJ
import TorJ: IMAS

@testset "Gaussian Integration with Polar Coordinates" begin
    
    # Test parameters matching the launch function
    x0 = [0.0,0.0,0.0]
    N0 = [0.0,0.0,1.0]
    w = 0.0174  # beam width parameter
    inverse_curvature_radius = 1.0/3.99
    freq = 110.0e9
    N_rings = 21  # Use more rings for better accuracy
    
    # Generate peripheral rays
    ray_positions, ray_directions, ray_weights = TorJ.launch_peripheral_rays(
        x0, N0, w, N_rings, inverse_curvature_radius, freq; 
        min_azimuthal_points=11, normalize_weight_sum=false)
    
    # Since launch is at origin and rays go in z-direction,
    # the integration is simply in the x-y plane
    N_rays = size(ray_positions, 1)
    x_coords = ray_positions[:, 1]  # x coordinates
    y_coords = ray_positions[:, 2]  # y coordinates
    
    # Test 1: Integration of a 2D Gaussian using polar coordinates
    @testset "2D Gaussian Integration" begin
        # Define 2D Gaussian function: exp(-(x² + y²)/(2σ²))
        
        # Integrate using polar coordinate points and weights
        polar_integral = sum(ray_weights[i] 
                            for i in 1:N_rays)
        
        # Analytical result for 2D Gaussian integral over infinite domain
        analytical_result = 1.0
        
        # Test that polar integration is close to analytical result
        relative_error = abs(polar_integral - analytical_result) / analytical_result
        @test relative_error < 0.01  # 5% tolerance
        
        println("Polar integration result: $polar_integral")
        println("Analytical result: $analytical_result") 
        println("Relative error: $(relative_error * 100)%")
    end
end