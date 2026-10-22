# ==============================================================================
# MOSTInversionBenchmark.jl
# Benchmark Script: C^inf Soft-Plus Solver vs. C^0 Unregularized max() Solver
# ==============================================================================

module MOSTInversionBenchmark

using BenchmarkTools
using ForwardDiff
using LinearAlgebra

# ------------------------------------------------------------------------------
# 1. Core Physics & Inversion Routines
# ------------------------------------------------------------------------------

struct MOSTParams{T<:AbstractFloat}
    Pr_t0::T
    β_m::T
    β_h::T
    Ri_c::T
    ϵ::T
    a1::T
    a2::T
    a3::T

    function MOSTParams{T}(; Pr_t0=1.0, β_m=4.7, β_h=4.7, Ri_c=0.19, ϵ=1e-5) where {T<:AbstractFloat}
        c1 = T(Pr_t0)
        c2 = T(β_h) - T(2.0) * T(β_m) * T(Pr_t0)
        c3 = T(3.0) * T(β_m)^2 * T(Pr_t0) - T(2.0) * T(β_m) * T(β_h)

        a1 = inv(c1)
        a2 = -c2 / (c1^3)
        a3 = (T(2.0) * c2^2 - c1 * c3) / (c1^5)

        new{T}(T(Pr_t0), T(β_m), T(β_h), T(Ri_c), T(ϵ), a1, a2, a3)
    end
end

MOSTParams() = MOSTParams{Float64}()

""" Unregularized C^0 max() Solver """
@inline function solve_zeta_max(Ri::T, p::MOSTParams{T}) where {T<:Real}
    abs_ri_eps = sqrt(Ri^2 + p.ϵ^2)
    Ri_pos = T(0.5) * (Ri + abs_ri_eps - p.ϵ)
    w_stable = T(0.5) * (T(1.0) + Ri / abs_ri_eps)

    p_ri = p.a1 + Ri * (p.a2 + Ri * p.a3)
    zeta_neutral = Ri * p_ri

    diff_c = T(1.0) - Ri_pos / p.Ri_c
    denom_stable = max(diff_c, p.ϵ)
    zeta_stable = (Ri_pos / p.Pr_t0) / denom_stable

    return (T(1.0) - w_stable) * zeta_neutral + w_stable * zeta_stable
end

""" C^inf Soft-Plus Smooth Solver """
@inline function solve_zeta_softplus(Ri::T, p::MOSTParams{T}) where {T<:Real}
    abs_ri_eps = sqrt(Ri^2 + p.ϵ^2)
    Ri_pos = T(0.5) * (Ri + abs_ri_eps - p.ϵ)
    w_stable = T(0.5) * (T(1.0) + Ri / abs_ri_eps)

    p_ri = p.a1 + Ri * (p.a2 + Ri * p.a3)
    zeta_neutral = Ri * p_ri

    diff_c = T(1.0) - Ri_pos / p.Ri_c
    denom_stable = T(0.5) * (diff_c + sqrt(diff_c^2 + p.ϵ^2))
    zeta_stable = (Ri_pos / p.Pr_t0) / denom_stable

    return (T(1.0) - w_stable) * zeta_neutral + w_stable * zeta_stable
end

# ------------------------------------------------------------------------------
# 2. Vector Evaluation Loops (Zero Allocation + SIMD)
# ------------------------------------------------------------------------------

function benchmark_vector_max!(out::Vector{T}, Ri_vec::Vector{T}, p::MOSTParams{T}) where {T}
    @inbounds @simd for i in eachindex(Ri_vec)
        out[i] = solve_zeta_max(Ri_vec[i], p)
    end
    return nothing
end

function benchmark_vector_softplus!(out::Vector{T}, Ri_vec::Vector{T}, p::MOSTParams{T}) where {T}
    @inbounds @simd for i in eachindex(Ri_vec)
        out[i] = solve_zeta_softplus(Ri_vec[i], p)
    end
    return nothing
end

# ------------------------------------------------------------------------------
# 3. Automatic Differentiation & Newton Inversion
# ------------------------------------------------------------------------------

function run_ad_jacobian_sweep(p::MOSTParams{Float64})
    Ri_threshold = p.Ri_c * (1.0 - p.ϵ)
    Ri_vals = range(Ri_threshold - 0.002, Ri_threshold + 0.002, length=500)

    d_max = [ForwardDiff.derivative(x -> solve_zeta_max(x, p), r) for r in Ri_vals]
    d_sp  = [ForwardDiff.derivative(x -> solve_zeta_softplus(x, p), r) for r in Ri_vals]

    kink_jump_max = maximum(abs.(diff(d_max)))
    smooth_jump_sp = maximum(abs.(diff(d_sp)))

    return kink_jump_max, smooth_jump_sp
end

function newton_inversion(solver_func, target_zeta::Float64, p::MOSTParams{Float64}; Ri_init=0.1899, max_iter=25, rel_tol=1e-6)
    Ri = Ri_init
    history = Tuple{Int, Float64, Float64, Float64}[]

    for iter in 1:max_iter
        z_val = solver_func(Ri, p)
        res = z_val - target_zeta
        j_val = ForwardDiff.derivative(x -> solver_func(x, p), Ri)

        push!(history, (iter - 1, Ri, res, j_val))

        # Relative residual check for well-conditioned Float64 termination
        if abs(res) / target_zeta < rel_tol
            return true, iter, history
        end

        j_safe = copysign(max(abs(j_val), 1e-12), j_val)
        step = res / j_safe
        Ri -= step
    end

    return false, max_iter, history
end

# ------------------------------------------------------------------------------
# 4. Main Driver Routine
# ------------------------------------------------------------------------------

function run_benchmark()
    p = MOSTParams()
    println("="^75)
    println("  MOSTInversionArchitectures: C^0 max() vs C^inf Soft-Plus Benchmark")
    println("="^75)

    N = 1_000_000
    Ri_vec = collect(range(-2.0, 0.25, length=N))
    out_max = zeros(N)
    out_sp  = zeros(N)

    println("Vector Size: N = $(N) points")
    println("\nEvaluating Vectorized Execution Throughput (@btime)...")
    print("  C^0 max() Kernel:       ")
    @btime benchmark_vector_max!($out_max, $Ri_vec, $p)

    print("  C^inf Soft-Plus Kernel: ")
    @btime benchmark_vector_softplus!($out_sp, $Ri_vec, $p)

    # ForwardDiff Sweep
    kink_max, smooth_sp = run_ad_jacobian_sweep(p)
    println("\nForwardDiff AD Jacobian Continuity near Ri ≈ Ri_c(1 - ε):")
    println("  C^0 max() Derivative Jump:       $(kink_max)  [C^0 KINK JUMP]")
    println("  C^inf Soft-Plus Derivative Jump: $(smooth_sp)  [SMOOTH C^inf GRADIENT]")

    # Newton-Raphson Inversion
    target_z = (p.Ri_c * (1.0 - p.ϵ) / p.Pr_t0) / p.ϵ
    conv_max, iters_max, hist_max = newton_inversion(solve_zeta_max, target_z, p)
    conv_sp, iters_sp, hist_sp   = newton_inversion(solve_zeta_softplus, target_z, p)

    println("\nNewton-Raphson Implicit Inversion near Critical Boundary:")
    println("  C^0 max() Status:      $(conv_max ? "CONVERGED" : "STALLED/FAILED") ($(iters_max) iters)")
    println("  C^inf Soft-Plus Status: $(conv_sp ? "CONVERGED" : "STALLED/FAILED") ($(iters_sp) iters)")

    println("\nFirst 5 Iteration Trajectory (C^inf Soft-Plus):")
    for (i, r, res, j) in hist_sp[1:min(5, length(hist_sp))]
        println("  Iter $i: Ri = $(round(r, digits=8)), Res = $(res), Jac = $(j)")
    end
    println("="^75)
end

end # module MOSTInversionBenchmark