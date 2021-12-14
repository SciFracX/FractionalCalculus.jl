import FractionalCalculus.FracDiffAlg

"""
Grünwald–Letnikov sense fractional derivative algorithms, please refer to [Grünwald–Letnikov derivative](https://en.wikipedia.org/wiki/Gr%C3%BCnwald%E2%80%93Letnikov_derivative) for more details
"""
abstract type GL <: FracDiffAlg end


"""
Grunwald Letnikov direct compute method to obtain fractional derivative, precision are guaranteed but cause more memory allocation and compilation time.
"""
struct GL_Direct <: GL end



"""
@book{oldham_spanier_1984,
title={The fractional calculus: Theory and applications of differentiation and integration to arbitrary order},
author={Oldham, Keith B. and Spanier, Jerome},
year={1984}} 
"""

"""
Using Grünwald–Letnikov multiplication-addition-multiplication-addition··· method to compute fractional derivative in Grünwald Letnikov sense.
"""
struct GL_Multiplicative_Additive <: GL end

"""
Using Lagrange three points interpolation method to compute fractional derivative in Grünwald Letnikov sense.
"""
struct GL_Lagrange_Three_Point_Interp <: GL end

"""
Using **Grünwald–Letnikov finite difference method** to compute Grünwald–Letnikov sense fractional derivative.
"""
struct GL_Finite_Difference <: GL end


################################################################
###                    Type defination done                  ###
################################################################


"""
# Grünwald–Letnikov sense fractional dervivative.

    fracdiff(f, α, start_point, end_point, GL_Direct())

### Example:

```julia-repl
julia> fracdiff(x->x^5, 0, 0.5, GL_Direct())
```

> Please note Grunwald-Letnikov sense fracdiff only support 0 < α < 1.

Please refer to [Grünwald–Letnikov derivative](https://en.wikipedia.org/wiki/Gr%C3%BCnwald%E2%80%93Letnikov_derivative) for more details.
"""
function fracdiff(f::Union{Function, Number}, α::Float64, start_point, end_point, ::GL_Direct)
    #checks(f, α, start_point, end_point)

    g(τ) = (f(end_point)-f(τ)) ./ (end_point - τ) .^ (1+α)
    result = f(end_point)/(gamma(1-α) .* (end_point - start_point) .^ α) .+ quadgk(g, start_point, end_point) .* α ./ gamma(1-α)
    return result
end

function fracdiff(f::Union{Number, Function}, α::Float64, start_point, end_point::AbstractArray, ::GL_Direct)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, start_point, value, GL_Direct()))
    end
    return ResultArray
end



#TODO: Use the improved alg!! This algorithm is not accurate
#This algorithm is not good, still more to do
"""
# Grünwald Letnikov sense derivative approximation

    fracdiff(f, α, end_point, h, GL_Multiplicative_Additive())

Grünwald–Letnikov multiplication-addition-multiplication-addition··· method to approximate fractional derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.007, GL_Multiplicative_Additive())
1.127403405642918
```

!!! danger "Inaccurate"
    The `GL_Multiplicative_Additive` method is not accruate, please use it at your own risk.
"""
function fracdiff(f::Union{Function, Number}, α, end_point, h, ::GL_Multiplicative_Additive)::Float64

    summation = 0
    n = Int64(floor(end_point/h))

    for i ∈ 0:n-1
        summation += gamma(i-α)/gamma(i+1)*f(end_point-i*h)
    end
    result = summation*end_point^(-α)*n^α/gamma(-α)

    return result
end

function fracdiff(f::Union{Number, Function}, α::Float64, end_point::AbstractArray, h, ::GL_Multiplicative_Additive)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, value, h, GL_Multiplicative_Additive()))
    end
    return ResultArray
end


#TODO: This algorithm is same with the above one, not accurate!!!
#This algorithm is not good, still more to do
"""
# Grünwald Letnikov sense derivative approximation

    fracdiff(f, α, end_point, h, GL_Lagrange_Three_Point_Interp())

Using Lagrange three poitns interpolation to approximate the fractional derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.007, GL_Lagrange_Three_Point_Interp())
1.1283963376811044
```

!!! danger "Inaccurate"
    The `GL_Multiplicative_Additive` method is not accruate, please use it at your own risk.
"""
function fracdiff(f::Union{Function, Number}, α::Float64, end_point, h, ::GL_Lagrange_Three_Point_Interp)::Float64
    #checks(f, α, 0, end_point)

    n = Int64(floor(end_point/h))
    summation=0

    for i ∈ 0:n-1
        summation += gamma(i-α)/gamma(i+1)*(f(end_point-i*h)+1/4*α*(f(end_point-(i-1)*h)-f(end_point-(i+1)*h))+1/8*α^2*(f(end_point-(i-1)*h)-2*f(end_point-i*h)+f(end_point-(i+1)*h)))
    end

    result = summation*end_point^(-α)*n^α/gamma(-α)

    return result
end

function fracdiff(f::Union{Number, Function}, α::Float64, end_point::AbstractArray, h, ::GL_Lagrange_Three_Point_Interp)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, value, h, GL_Lagrange_Three_Point_Interp()))
    end
    return ResultArray
end


"""
# Grünwald Letnikov sense derivative approximation

    fracdiff(f::Union{Function, Number}, α::AbstractArray, end_point, h, ::GL_Finite_Difference)::Vector

Use finite difference method to obtain Grünwald Letnikov sense fractional derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.01, GL_Finite_Difference())
1.1269695801851276
```

!!! danger "Inaccurate"
    `GL_Finite_Difference` method is not accruate, please use it at your own risk.
"""
function fracdiff(f::Union{Number, Function}, α::Float64, end_point::Real, h, ::GL_Finite_Difference)::Float64

    n = Int64(floor(end_point/h))
    result=0

    @fastmath @simd for i ∈ 0:n
        result += (-1)^i/(gamma(i+1)*gamma(α-i+1))*f(end_point-i*h)
    end

    result1=result*h^(-α)*gamma(α+1)
    return result1
end

function fracdiff(f::Union{Function, Number}, α::AbstractArray, end_point, h, ::GL_Finite_Difference)::Vector
    ResultArray = Float64[]

    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracint(f, α, value, h, GL_Finite_Difference()))
    end
    return ResultArray    
end