"""

Grünwald–Letnikov sense fractional derivative algorithms, please refer to [Grünwald–Letnikov derivative](https://en.wikipedia.org/wiki/Gr%C3%BCnwald%E2%80%93Letnikov_derivative) for more details
"""
abstract type GL <: FracDiffAlg end


"""
# Grünwald–Letnikov sense fractional dervivative.

    fracdiff(f, α, start_point, end_point, GLDirect())

### Example:

```julia-repl
julia> fracdiff(x->x^5, 0, 0.5, GLDirect())
```

!!! info "Scope"
    Please note Grunwald-Letnikov sense fracdiff only support ``0 < \\alpha < 1``.

Please refer to [Grünwald–Letnikov derivative](https://en.wikipedia.org/wiki/Gr%C3%BCnwald%E2%80%93Letnikov_derivative) for more details.

Grunwald Letnikov direct compute method to obtain fractional derivative, precision are guaranteed but cause more memory allocation and compilation time.
"""
struct GLDirect <: GL end




"""
# Grünwald Letnikov sense Multiplicative and Addtive approximation

    fracdiff(f, α, end_point, h, GLMultiplicativeAdditive())

Grünwald–Letnikov multiplication-addition-multiplication-addition··· method to approximate fractional derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.007, GLMultiplicativeAdditive())
1.127403405642918
```

```tex
@book{oldham_spanier_1984,
title={The fractional calculus: Theory and applications of differentiation and integration to arbitrary order},
author={Oldham, Keith B. and Spanier, Jerome},
year={1984}} 
```
"""
struct GLMultiplicativeAdditive <: GL end



"""
# Grünwald Letnikov sense three point interpolation

    fracdiff(f, α, end_point, h, GLLagrangeThreePointInterp())

Using Lagrange three poitns interpolation to approximate the fractional derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.006, GLLagrangeThreePointInterp())
1.1261297605404632
```
"""
struct GLLagrangeThreePointInterp <: GL end

"""
# Grünwald Letnikov sense finite difference approximation

    fracdiff(f::Union{Function, Number}, α::AbstractArray, end_point, h, ::GLFiniteDifference)::Vector

Use finite difference method to obtain Grünwald Letnikov sense fractional derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.01, GLFiniteDifference())
1.1269695801851276
```
"""
struct GLFiniteDifference <: GL end

"""
# Grünwald Letnikov sense derivative approximation

    fracdiff(f, α, point, p, GLHighPrecision())

Use the high precision algorithms to compute the Grunwald letnikov fractional derivative.

The **p** here is the grade of precision.

!!! note
    The value interval passing in the function should be a array!
"""
struct GLHighPrecision <: GL
    p::Int
end

GLHighPrecision(; p::Int=2) = GLHighPrecision(p)


################################################################
###                    Type definition done                  ###
################################################################


#=
Grunwald Letnikov direct method
=#
function fracdiff(f::FunctionAndNumber,
                  α::T,
                  start_point::Real,
                  end_point::Real,
                  ::GLDirect) where {T <: Real}
    #checks(f, α, start_point, end_point)
    typeof(f) <: Number ? (end_point == 0 ? (return 0) : (return f/sqrt(pi*end_point))) : nothing
    end_point == 0 ? (return 0) : nothing
    g(τ) = (f(end_point)-f(τ)) ./ (end_point - τ) .^ (1+α)
    result = f(end_point)/(gamma(1-α) .* (end_point - start_point) .^ α) .+ (quadgk(g, start_point, end_point) .* α ./ gamma(1-α))[1]
    return result
end

fracdiff(f::FunctionAndNumber, α::Float64, end_point, ::GLDirect) = fracdiff(f::FunctionAndNumber, α::Float64, 0, end_point, GLDirect())


#TODO: Use the improved alg!! This algorithm is not accurate
#This algorithm is not so good, still more to do
function fracdiff(f::FunctionAndNumber, α, end_point, h::Real, ::GLMultiplicativeAdditive)
    typeof(f) <: Number ? (end_point == 0 ? 0 : f/sqrt(pi*end_point)) : nothing
    summation = zero(Float64)
    n = round(Int, end_point/h)

    for i ∈ 0:n-1
        summation += gamma(i-α)/gamma(i+1)*f(end_point-i*h)
    end
    result = summation*end_point^(-α)*n^α/gamma(-α)

    return result
end
#=
function fracdiff(f::Union{Number, Function}, α::Float64, end_point::AbstractArray, h, ::GLMultiplicativeAdditive)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, value, h, GLMultiplicativeAdditive()))
    end
    return ResultArray
end
=#

#TODO: This algorithm is same with the above one, not accurate!!!
#This algorithm is not so good, still more to do
function fracdiff(f::FunctionAndNumber,
                  α::Real,
                  end_point::Real,
                  h::Real,
                  ::GLLagrangeThreePointInterp)
    #checks(f, α, 0, end_point)
    typeof(f) <: Number ? (end_point == 0 ? 0 : f/sqrt(pi*end_point)) : nothing
    n = round(Int, end_point/h)
    summation = zero(Float64)

    for i ∈ 0:n-1
        summation += gamma(i-α)/gamma(i+1)*(f(end_point-i*h)+1/4*α*(f(end_point-(i-1)*h)-f(end_point-(i+1)*h))+1/8*α^2*(f(end_point-(i-1)*h)-2*f(end_point-i*h)+f(end_point-(i+1)*h)))
    end

    result = summation*end_point^(-α)*n^α/gamma(-α)
    return result
end
#=
function fracdiff(f::Union{Number, Function}, α::Float64, end_point::AbstractArray, h, ::GLLagrangeThreePointInterp)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, value, h, GLLagrangeThreePointInterp()))
    end
    return ResultArray
end
=#

function fracdiff(f::FunctionAndNumber,
                  α::Real,
                  end_point::Real,
                  h::Real,
                  ::GLFiniteDifference)
    typeof(f) <: Number ? (end_point == 0 ? 0 : f/sqrt(pi*end_point)) : nothing
    n = round(Int, end_point/h)
    result = zero(Float64)

    @fastmath @simd for i = 0:n
        result += (-1)^i/(gamma(i+1)*gamma(α-i+1))*f(end_point-i*h)
    end

    result1 = result*h^(-α)*gamma(α+1)
    return result1
end
#=
function fracdiff(f::Union{Function, Number}, α::AbstractArray, end_point, h, ::GLFiniteDifference)::Vector
    ResultArray = Float64[]

    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracint(f, α, value, h, GLFiniteDifference()))
    end
    return ResultArray    
end
=#



function fracdiff(f::Union{Function, Number, Vector},
                  α::T,
                  t,
                  alg::GLHighPrecision) where {T <: Real}
    @unpack p = alg
    if isa(f, Function)
        y=f.(t)
    elseif isa(f, Vector)
        y=f[:]
    end
    h = t[2] - t[1]
    
    t = t[:]
    n = length(t)
    u = 0
    du = 0
    r = collect(0:p)*h
    R = reverse(Vandermonde(r), dims=2)
    c = inv(R)*y[1:p+1]
    for i = 1:p+1
        u = u.+c[i]*t.^(i-1)
        du = du.+c[i]*t.^(i-1-α)*gamma(i)/gamma(i-α)
    end

    v = y.-u
    g = genfun(p)
    w = getvec(α, n, g)
    dv = zeros(n)
    for i=1:n
        dv[i] = w[1:i]'*v[i:-1:1]/h^α
    end

    dy = dv .+ du
    if abs(y[1])<1e-10
        dy[1]=0
    end
    return dy
end

"""
P-th precision polynomial generate function

```math
g_p(z)=\\sum_{k=1}^p \\frac{1}{k}(1-z)^k
```
"""
function genfun(p)
    a = collect(1:p+1)
    A = Vandermonde(a)'
    return (1 .-a')*inv(A')
end

function getvec(α, n, g)
    p = length(g)-1
    b = 1 + α
    g0 = g[1]
    w = Float64[]
    push!(w, g[1]^α)

    for m = 2:p
        M = m-1
        dA = b/M
        temp = (-(g[2:m] .*collect((1-dA):-dA:(1-b))))' *w[M:-1:1]/g0
        push!(w, temp)
    end

    for k = p+1:n
        M = k-1
        dA = b/M
        temp = (-(g[2:(p+1)] .*collect((1-dA):-dA:(1-p*dA))))' *w[M:-1:(k-p)]/g0
        push!(w, temp)
    end

    return w
end