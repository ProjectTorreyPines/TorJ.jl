import Enzyme: autodiff, Duplicated, Const, Reverse

function gradΛ_enzyme!(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, plasma::Plasma, ω::Real, mode::Integer)
    # Struggling to get rid of this allocation.
    # It should be easy to use `du` as a temporary array but 
    # in practice it always leads to errors
    ∂Λ_∂u = zeros(Float64, 6)
    autodiff(Reverse, dispersion_relation,
             Duplicated(u, ∂Λ_∂u),      # Compute ∂Λ/∂x
             Const(plasma),
             Const(ω),
             Const(mode))
    # Normalize with | ∂Λ_∂N |
    norm_∂Λ_∂N = LinearAlgebra.norm(∂Λ_∂u[4:6])
    # # Swap du[1:3] with du[4:6]
    # # dx/ds = ∂Λ_∂N and dN/ds = ∂Λ_∂x
    du[1:3] .= ∂Λ_∂u[4:6] /norm_∂Λ_∂N
    du[4:6] .= -∂Λ_∂u[1:3] /norm_∂Λ_∂N
end