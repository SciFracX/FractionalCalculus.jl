"""
Hadamard sense fractional derivative algorithms.
"""
abstract type HadamardDiff <: FracDiffAlg end


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
# Hadamard sense left rectangular approximating algorithm

    fracdiff(f, α, a, x, h, HadamardLRect())

Using finite part integral of **left rectangular formula** to compute the fractional hadamard derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 0, 1, 0.001, HadamardLRect())
0.9738507879228357
```
"""
struct HadamardLRect <: HadamardDiff end

"""
# Hadamard sense right rectangular approximating algorithm

    fracdiff(f, α, a, x, h, HadamardTrap())

Using finite part integral of **right rectangular formula** to compute the fractional hadamard derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 0, 1, 0.001, HadamardRRect())
0.9738507879228337
```
"""
struct HadamardRRect <: HadamardDiff end

"""
# Hadamard sense trapezoidal approximating algorithm

    fracdiff(f, α, a, x, h, HadamardTrap())

Using finite part integral of **trapezoidal formula** to compute the fractional hadamard derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 0, 1, 0.001, HadamardTrap())
0.9738507879228328
```
"""
struct HadamardTrap <: HadamardDiff end


################################################################
###                    Type definition done                  ###
################################################################

# FIXME: Need to verify the derivative value in the interval [x₀, x]
#=
Hadamard Left Rectangular computing algorithms
=#
function fracdiff(f::FunctionAndNumber,
                  α::Float64,
                  x₀,
                  x::Real,
                  h::Float64,
                  ::HadamardLRect)
    typeof(f) <: Number ? (x == 0 ? (return 0) : (return f/sqrt(pi*x))) : nothing
    x == 0 ? (return 0) : nothing
    N = round(Int, (x-x₀)/h)
    result = zero(Float64)

    @fastmath @inbounds @simd for i ∈ 0:N-1
        result += LRectCoeff(i, N, h, α, x₀)*f(x₀+i*h)
    end

    return 1/gamma(1-α)*result
end
function LRectCoeff(i::Int64, n::Int64, h::Float64, α::Float64, x₀)
    if 0 ≤ i ≤ n-2
        return ((log((x₀+n*h)/(x₀+i*h)))^(-α) - (log((x₀+n*h)/(x₀+(i+1)*h)))^(-α))
    elseif i == n-1
        return (log((x₀+n*h)/(x₀+i*h)))^(-α)
    end
end

#=
Hadamard Right Rectangular computing algorithm
=#
function fracdiff(f::FunctionAndNumber,
                  α::Float64,
                  x₀::Real,
                  x::Real,
                  h::Float64,
                  ::HadamardRRect)
    typeof(f) <: Number ? (x == 0 ? 0 : f/sqrt(pi*x)) : nothing
    N = round(Int, (x-x₀)/h)
    result = zero(Float64)

    @fastmath @inbounds @simd for i ∈ 0:N-1
        result += LRectCoeff(i, N, h, α, x₀)*f(x₀+(i+1)*h)
    end

    return 1/gamma(1-α)*result
end
function RRectCoeff(i::Int64, n::Int64, h::Float64, α::Float64, x₀::Real)
    if i == 0
        return 1/2*LRectCoeff(i, n, h, α, x₀)
    elseif 1 ≤ i ≤ n-1
        return 1/2*(LRectCoeff(i-1, n, h, α, x₀) + LRectCoeff(i, n, h, α, x₀))
    elseif i == n
        return 1/2*LRectCoeff(i-1, n, h, α, x₀)
    end
end

#=
Hadamard trapezoidal computing algorithm
=#
function fracdiff(f::FunctionAndNumber,
                  α::Float64,
                  x₀::Real,
                  x::Real,
                  h::Float64,
                  ::HadamardTrap)
    typeof(f) <: Number ? (x == 0 ? 0 : f/sqrt(pi*x)) : nothing
    N = round(Int, (x-x₀)/h)

    result = zero(Float64)

    @fastmath @inbounds @simd for i ∈ 0:N
        result += RRectCoeff(i, N, h, α, x₀)*f(x₀+i*h)
    end

    return 1/gamma(1-α)*result
end

function fracdiff(f, α, x, h, ::HadamardTrap) = fracdiff(f, α, 0, x, h, HadamardTrap())