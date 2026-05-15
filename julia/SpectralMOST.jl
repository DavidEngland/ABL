module SpectralMOST

export stability_profile, integral_correction

"""
    stability_profile(ζ, λ; N=20)

Computes the non-dimensional gradient function ϕ(ζ) = 1 + Σ α_n (γ ζ)^n
using the recursive spectral expansion.
λ = 1/4 for Momentum (m), λ = 1/2 for Heat (h).
"""
function stability_profile(ζ::T, λ::T; N::Int=20, γ=15.0) where T<:AbstractFloat
    # Start with n=0 (α_0 = 1)
    # We skip the '1' in the sum and add it at the end
    α_n = 1.0
    total = 0.0
    γζ = γ * ζ

    for n in 1:N
        # Recursive step: α_n = α_{n-1} * (n - (1-λ)) / n
        α_n *= (n - (1 - λ)) / n
        total += α_n * (γζ^n)
    end

    return 1.0 + total
end

"""
    integral_correction(ζ, λ; N=20)

Computes the integrated stability correction ψ(ζ) = Σ β_n (γ ζ)^n
where β_n = α_n / n. This reweights the spectrum toward low modes.
"""
function integral_correction(ζ::T, λ::T; N::Int=20, γ=15.0) where T<:AbstractFloat
    α_n = 1.0
    total = 0.0
    γζ = γ * ζ

    for n in 1:N
        α_n *= (n - (1 - λ)) / n
        # The Integral Damping: β_n = α_n / n
        β_n = α_n / n
        total += β_n * (γζ^n)
    end

    return total
end

end # module

# --- Usage Example ---
using .SpectralMOST

ζ = 0.1
ψ_m = integral_correction(ζ, 0.25) # Momentum (λ=1/4)
ψ_h = integral_correction(ζ, 0.50) # Heat (λ=1/2)

println("At ζ=$ζ:")
println("ψ_m (Spectral): ", ψ_m)
println("ψ_h (Spectral): ", ψ_h)