const fast = 3
const i_max = 3

function warmdisp(X::Real, Y::Real, N_par::Real, mu::Real, imod:: Integer, fast:: Integer, lrm:: Integer, Npr::Real)
    return 0.0
end

function absorption(X::Real, Y::Real, N_r::Real, theta::Real, te:: Real, Real:: v_g_perp, imod:: Integer)
    N_par = Nr * cos(theta)
    Npr = sqrt(N_r^2 - N_par^2)
    mu = constants.m_e * constants.c^2/(te * e0)
    nharm = larmornumber(Y, Nll, mu)
    lrm = min(i_max, nharm)
    N_perp_cmplx = warmdisp(X, Y, N_par, mu, imod, fast, lrm, Npr)
    return 2.e0 * imag(N_perp_cmplx^2) * omega / c0 * v_g_perp
end