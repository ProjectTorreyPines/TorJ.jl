const ntv = 501
const tmax = 5.0
const dt = 2.0 * tmax / (ntv - 1)
const i_max = 5 

# Module-level variables (equivalent to Fortran save variables)
_int_weights = Float64[]
_int_absz = Float64[]

_ttv = Vector{Float64}(undef, ntv)
_extdtv = Vector{Float64}(undef, ntv)

const constants = (
    μ_0=1.25663706212e-6,
    c=2.99792458e8,
    ϵ_0=8.8541878128e-12,
    k_B=1.380649e-23,
    e=1.602176634e-19,
    m_e=9.1093837015e-31,
    m_p=1.67262192369e-27,
    m_n=1.67492749804e-27,
    atm=101325.0,
    m_u=1.6605390666e-27,
    avog=6.02214076e23,
    π_sqrt=1.7724538509055160272981674833411
)