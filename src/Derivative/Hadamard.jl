import FractionalCalculus.FracDiffAlg

"""
Hadamard sense fractional derivative algorithms.
"""
abstract type Hadamard <: FracDiffAlg end


"""
@article{article,
author = {Yin, Chuntao and Li, Changpin and Bi, Qinsheng},
year = {2018},
month = {01},
pages = {},
title = {Approximating Hadamard Derivative and Fractional Differential Equation via the Finite Part Integral},
journal = {SSRN Electronic Journal},
doi = {10.2139/ssrn.3281675}
}
"""

"""
Using finite part integral of **left rectangular formula** to compute the fractional hadamard derivative.
"""
struct Hadamard_LRect <: Hadamard end

"""
Using finite part integral of **right rectangular formula** to compute the fractional hadamard derivative.
"""
struct Hadamard_RRect <: Hadamard end

"""
Using finite part integral of **trapezoidal formula** to compute the fractional hadamard derivative.
"""
struct Hadamard_Trap <: Hadamard end


################################################################
###                    Type defination done                  ###
################################################################



"""
    fracdiff(f, α, a, x, h, Hadamard_LRect())

Compute Hadamard sense using left rectangular formula.
"""
function fracdiff(f, α, x₀, x, h, ::Hadamard_LRect)
    N = Int64((x-x₀)/h)
    result = 0

    for i ∈ 0:N-1
        result += LRectCoeff(i, N, h, α, x₀)*f(x₀+i*h)
    end

    return result
end
function LRectCoeff(i, n, h, α, x₀)
    if 0 ≤ i ≤ n-2
        return 1/gamma(1-α)*((log((x₀+n*h)/(x₀+i*h)))^(-α) - (log((x₀+n*h)/(x₀+(i+1)*h)))^(-α))
    elseif i == n-1
        return 1/gamma(1-α)*(log((x₀+n*h)/(x₀+i*h)))^(-α)
    end
end

"""
    fracdiff(f, α, a, x, h, Hadamard_Trap())

Compute Hadamard sense using right rectangular formula.
"""
function fracdiff(f, α, x₀, x, h, ::Hadamard_RRect)
    N = Int64((x-x₀)/h)
    result = 0

    for i ∈ 0:N-1
        result += LRectCoeff(i, N, h, α, x₀)*f(x₀+(i+1)*h)
    end

    return result
end
function RRectCoeff(i, n, h, α, x₀)
    if i == 0
        return 1/2*LRectCoeff(i, n, h, α, x₀)
    elseif 1 ≤ i ≤ n-1
        return 1/2*(LRectCoeff(i-1, n, h, α, x₀)+LRectCoeff(i, n, h, α, x₀))
    elseif i == n
        return 1/2*LRectCoeff(i-1, n, h, α, x₀)
    end
end

"""
    fracdiff(f, α, a, x, h, Hadamard_Trap())

Compute Hadamard sense using Trapezoidal formula.
"""
function fracdiff(f, α, x₀, x, h, ::Hadamard_Trap)
    N=Int64((x-x₀)/h)

    result=0

    for i ∈ 0:N
        result += RRectCoeff(i, N, h, α, x₀)*f(x₀+i*h)
    end

    return result
end