function abs_Al_init(N_absz::Int)
    # Get Gauss-Legendre quadrature points and weights on [-1, 1]
    global _int_weights, _int_absz
    _int_absz, _int_weights = gausslegendre(N_absz)
    _int_absz = Vector{Float64}(_int_absz)
    _int_weights = Vector{Float64}(_int_weights)
end


function abs_Al_N_with_pol_vec(X::Float64, Y::Float64, cos_theta::Float64, sin_theta::Float64, mode::Int)
    # mode==1 --> X-mode, mode=-1 --> O-mode
    # Returns N and polarization vector e
    e = zeros(ComplexF64, 3)
    
    if X >= 1.0
        return 0.0, e
    end
    
    rho = Y^2 * sin_theta^4 + 4.0 * (1.0 - X)^2 * cos_theta^2
    if rho < 0.0
        return 0.0, e
    end
    
    rho = sqrt(rho)
    f = (2.0 * (1.0 - X)) / (2.0 * (1.0 - X) - Y^2 * sin_theta^2 - Float64(mode) * Y * rho)
    N = 1.0 - X * f
    
    if N < 0.0
        return 0.0, e
    end
    
    N = sqrt(N)
    
    if cos_theta^2 < 1e-5 || 1.0 - sin_theta^2 < 1e-5
        # Quasi-perpendicular
        if mode > 0  # X-mode
            e[2] = complex(0.0, sqrt(1.0 / N))
            e[1] = complex(0.0, 1.0 / Y * (1.0 - (1.0 - Y^2) * f)) * e[2]
            # e_z zero for quasi-perpendicular X-mode
        else
            e[3] = sqrt(1.0 / N)
        end
    else
        a_sq = sin_theta^2 * (1.0 + (((1.0 - X) * N^2 * cos_theta^2) / 
            (1.0 - X - N^2 * sin_theta^2)^2) * 
            1.0 / Y^2 * (1.0 - (1.0 - Y^2) * f)^2)^2
        b_sq = cos_theta^2 * (1.0 + ((1.0 - X) / 
            (1.0 - X - N^2 * sin_theta^2)) * 
            1.0 / Y^2 * (1.0 - (1.0 - Y^2) * f)^2)^2
        
        if mode > 0
            # X-mode - e_y positive, purely imaginary
            e[2] = complex(0.0, sqrt(1.0 / (N * sqrt(a_sq + b_sq))))
        else
            # O-mode - e_y negative, purely imaginary
            e[2] = complex(0.0, -sqrt(1.0 / (N * sqrt(a_sq + b_sq))))
        end
        
        e[1] = complex(0.0, 1.0 / Y * (1.0 - (1.0 - Y^2) * f)) * e[2]
        e[3] = -complex((N^2 * sin_theta * cos_theta) / (1.0 - X - N^2 * sin_theta^2), 0.0) * e[1]
    end
    
    return N, e
end

function get_upper_limit_tau(mu::Float64, Y::Float64, N_par::Float64, omega::Float64, ds::Float64)
    omega_bar = 1.0 / Y
    m_max = 3
    f_sum = 0.0
    
    for m in 2:m_max
        m_omega_bar = Float64(m) / omega_bar
        if m_omega_bar^2 + N_par^2 < 1.0
            continue
        end
        
        if m == 2
            t1 = (-8 * exp((mu * (-1 + N_par^2 + m_omega_bar - N_par * sqrt(-1 + N_par^2 + m_omega_bar^2))) / 
                (-1 + N_par^2)) * (-3 - N_par^4 * (3 + mu^2) + 
                3 * N_par * mu * sqrt(-1 + N_par^2 + m_omega_bar^2) - 
                3 * N_par^3 * mu * sqrt(-1 + N_par^2 + m_omega_bar^2) + 
                N_par^2 * (6 - mu^2 * (-1 + m_omega_bar^2)))) / (N_par^5 * mu^5)
            
            t2 = (-8 * exp((mu * (-1 + N_par^2 + m_omega_bar + N_par * sqrt(-1 + N_par^2 + m_omega_bar^2))) / 
                (-1 + N_par^2)) * (3 + N_par^4 * (3 + mu^2) + 
                3 * N_par * mu * sqrt(-1 + N_par^2 + m_omega_bar^2) - 
                3 * N_par^3 * mu * sqrt(-1 + N_par^2 + m_omega_bar^2) + 
                N_par^2 * (-6 + mu^2 * (-1 + m_omega_bar^2)))) / (N_par^5 * mu^5)
        else
            t1 = (48 * exp((mu * (-1 + N_par^2 + m_omega_bar - N_par * sqrt(-1 + N_par^2 + m_omega_bar^2))) / 
                (-1 + N_par^2)) * (-15 + 3 * N_par^6 * (5 + 2 * mu^2) + 
                15 * N_par * mu * sqrt(-1 + N_par^2 + m_omega_bar^2) + 
                N_par^5 * mu * (15 + mu^2) * sqrt(-1 + N_par^2 + m_omega_bar^2) + 
                3 * N_par^4 * (-15 + 2 * mu^2 * (-2 + m_omega_bar^2)) + 
                N_par^2 * (45 - 6 * mu^2 * (-1 + m_omega_bar^2)) + 
                N_par^3 * mu * sqrt(-1 + N_par^2 + m_omega_bar^2) * (-30 + mu^2 * (-1 + m_omega_bar^2)))) / 
                (N_par^7 * mu^7)
            
            t2 = (48 * exp((mu * (-1 + N_par^2 + m_omega_bar + N_par * sqrt(-1 + N_par^2 + m_omega_bar^2))) / 
                (-1 + N_par^2)) * (15 - 3 * N_par^6 * (5 + 2 * mu^2) + 
                15 * N_par * mu * sqrt(-1 + N_par^2 + m_omega_bar^2) + 
                N_par^5 * mu * (15 + mu^2) * sqrt(-1 + N_par^2 + m_omega_bar^2) + 
                N_par^4 * (45 - 6 * mu^2 * (-2 + m_omega_bar^2)) + 
                N_par^3 * mu * sqrt(-1 + N_par^2 + m_omega_bar^2) * (-30 + mu^2 * (-1 + m_omega_bar^2)) + 
                N_par^2 * (-45 + 6 * mu^2 * (-1 + m_omega_bar^2)))) / (N_par^7 * mu^7)
        end
        
        f = (t1 + t2) * mu * (N_par^2 + 1.0) * m^2
        
        if m == 3
            f = f * 20.25 * (1.0 - N_par^2)^2 * 0.5^6
        else
            f = f * 4.0 * (1.0 - N_par^2) * 0.5^4
        end
        
        f_sum = f_sum + f
    end
    
    if f_sum == 0.0
        return 0.0
    end
    
    norm = 1.0 / (1.0 + 105.0 / (128.0 * mu^2) + 15.0 / (8.0 * mu))
    result = f_sum * norm * (sqrt(mu / (2 * π)))^3 * ds
    result = result * (constants.e * Y)^2 * π * 8.0 * omega
    
    return result
end



function abs_Al_pol_fact(t::Vector{Float64}, omega_bar::Float64, m_0::Float64, 
                         N_par::Float64, N_perp::Float64, e::Vector{ComplexF64}, m::Int)
    
    x_m = N_perp * omega_bar * sqrt((Float64(m) / m_0)^2 - 1.0)
    N_eff = (N_perp * N_par) / (1.0 - N_par^2)
    
    pol_vect = copy(e)
    
    Axz = e[1] + N_eff * e[3]
    Axz_sq = abs(Axz)^2  # mixed terms do not vanish!
    Re_Axz_ey = real(complex(0.0, 1.0) * Axz * conj(pol_vect[2]))
    Re_Axz_ez = real(Axz * conj(pol_vect[3]))
    Re_ey_ez = real(complex(0.0, 1.0) * conj(pol_vect[2]) * pol_vect[3])
    ey_sq = abs(pol_vect[2])^2
    ez_sq = abs(pol_vect[3])^2
    
    bessel_arg = x_m * sqrt.(1.0 .- t.^2)
    
    # Note: You'll need to implement or import Bessel functions
    # Using SpecialFunctions.jl package: besselj(n, x)
    bessel_n_l = [besselj(m - 1, arg) for arg in bessel_arg]
    bessel_n = [besselj(m, arg) for arg in bessel_arg]
    bessel_n_2 = bessel_n.^2
    bessel_n_u = [besselj(m + 1, arg) for arg in bessel_arg]
    
    abs_Al_bessel_sqr_deriv = bessel_arg ./ x_m .* bessel_n .* (bessel_n_l .- bessel_n_u)
    
    pol_fact = (Axz_sq + ey_sq) * bessel_n_2
    pol_fact .+= Re_Axz_ey * x_m / Float64(m) * abs_Al_bessel_sqr_deriv
    pol_fact .-= (bessel_arg ./ Float64(m)).^2 .* ey_sq .* bessel_n_l .* bessel_n_u
    pol_fact .+= (x_m / (Float64(m) * sqrt(1.0 - N_par^2)))^2 .* ez_sq .* t.^2 .* bessel_n_2
    pol_fact .+= x_m / (Float64(m) * sqrt(1.0 - N_par^2)) .* 2.0 .* Re_Axz_ez .* t .* bessel_n_2
    pol_fact .+= x_m / (Float64(m) * sqrt(1.0 - N_par^2)) .* Re_ey_ez .* t .* x_m / Float64(m) .* abs_Al_bessel_sqr_deriv
    pol_fact .*= (Float64(m) / (N_perp * omega_bar))^2
    
    return pol_fact
end

function abs_Al_integral_nume_fast(mu::Float64, omega_bar::Float64, m_0::Float64, 
                                   N_par::Float64, N_perp::Float64, e::Vector{ComplexF64}, m::Int)
    
    if length(_int_weights) == 0
        throw(ErrorException("The weights and abscissae for the absorption were never initialized. Call `abs_Al_init` before using the absorption."))
    end
    u_par = 1.0 / sqrt(1.0 - N_par^2) .* (Float64(m) / m_0 * N_par .+ 
            sqrt((Float64(m) / m_0)^2 - 1.0) .* _int_absz)
    u_perp_sq = ((Float64(m) / m_0)^2 - 1.0) .* (1.0 .- _int_absz.^2)
    gamma = sqrt.(1.0 .+ u_par.^2 .+ u_perp_sq)
    pol_fact = abs_Al_pol_fact(_int_absz, omega_bar, m_0, N_par, N_perp, e, m)
    
    c_abs_int = _int_weights .* pol_fact .* (-mu) .* exp.(mu .* (1.0 .- gamma))
    c_abs = sum(c_abs_int)
    
    a = 1.0 / (1.0 + 105.0 / (128.0 * mu^2) + 15.0 / (8.0 * mu))
    c_abs = c_abs * a * (sqrt(mu / (2 * π)))^3
    
    return c_abs
end

function abs_Albajar_fast(omega::Float64, X::Float64, Y::Float64, N_abs::Float64, 
                          N_par::Float64, Te::Float64, mode::Int)
    
    if Te < 20.0
        return 0.0  # very low Te => absorption can be ignored
    end
    
    mu = constants.m_e * constants.c^2 / (constants.e * Te)
    max_harmonic = 3  # consider only second and third harmonic
    omega_bar = 1.0 / Y
    c_abs = 0.0
    cos_theta = N_par/N_abs
    sin_theta = sin(acos(cos_theta))
    N_perp = sqrt(N_abs^2 - N_par^2)
    
    N_abs_test, e_pol = abs_Al_N_with_pol_vec(X, Y, cos_theta, sin_theta, mode)
    if isnan(N_abs_test) || N_abs_test <= 0.0 || N_abs_test > 1.0
        return 0.0
    end
    
    m_0 = sqrt(1.0 - N_par^2) * omega_bar
    
    for m_sum in 2:max_harmonic  # First harmonic needs to be treated separately (for now ignored)
        if Float64(m_sum) < m_0
            continue
        end
        c_abs_m = abs_Al_integral_nume_fast(mu, omega_bar, m_0, N_par, N_perp, e_pol, m_sum)
        c_abs += sqrt((Float64(m_sum) / m_0)^2 - 1.0) * c_abs_m
    end
    
    c_abs = -(c_abs * 2.0 * π^2 / m_0)  # Splitting this is just for overview
    
    c_abs = c_abs * X * omega / (Y * constants.c)  # revert the normalization
    
    return c_abs
end

function α_approx(x::AbstractVector{<:Real}, N::AbstractVector{<:Real}, plasma::Plasma, 
                  omega:: Real, mode::Integer)
    N_abs = LinearAlgebra.norm(N)
    X, Y, N_par, b = eval_plasma(plasma, x, N, omega)
    Te = T_e(plasma, x)
    α = abs_Albajar_fast(omega, X, Y, N_abs, N_par, Te, mode)
    return α
end