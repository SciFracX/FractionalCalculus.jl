"""
# Riemann-Liouville sense fractional integral algorithms

As the general fractional order integral type, RLInt is the parent type of all Riemann Liouville sense numerical methods.
"""
abstract type RLInt <: FracIntAlg end

"""
# Riemann Liouville sense direct compute

    fracint(f::Function, α, start_point, end_point, h, RLDirect())

Riemann_Liouville fractional integral using complex step differentiation.

```math
(J^αf)(t)=\\frac{1}{\\Gamma(α)} \\int_0^t(t-τ)^{α-1}f(τ)dτ
```

By using [QuadGK](https://github.com/JuliaMath/QuadGK.jl) calculate the integration and obtain the value.

### Example: 

```julia-repl
julia> fracint(x->x^5, 0.5, 0, 2.5, 1e-8, RLDirect())
(64.36234189727955, 7.742994054849976e-7)
```

Returns a tuple (1.1639316474512205, 1.0183453796725215e-8), which contains the value of this derivative is 1.1639316474512205, and the error estimate is 1.0183453796725215e-8
"""
struct RLDirect <: RLInt end
#=
"""
# Riemann Liouville sense fractional integral with first diff known.

    fracint(f, fd, α, start_point, end_point, RLDirect_First_Diff_Known())

### Example

```julia-repl
julia> fracint(x->x^5, x->5*x^4, 0.5, 0, 2.5, RLDirect_First_Diff_Known())
```

With first derivative known, we can use the Riemann Liouville sense to obtain the fractional integral more effcient.
"""
struct RLDirect_First_Diff_Known <: RLInt end
=#
"""
@article{LI20113352,
title = {Numerical approaches to fractional calculus and fractional ordinary differential equation},
author = {Changpin Li and An Chen and Junjie Ye},
}
"""

"""
# Riemann Liouville sense fractional integral using piecewise interpolation.
    
    fracint(f, α, end_point, h, RLPiecewise())
    
### Example
    
```julia-repl
julia> fracint(x->x^5, 0.5, 2.5, 0.0001, RLPiecewise())
64.36234206345209
```
    
By deploying Piecewise interpolation to approximate the original function, with small step size, this method is fast and take little memory allocation.

Using piecewise linear interpolation:

```math
    y_n(t)=\\frac{t-t_{i+1}}{t_i-t_{i+1}}y(t_i)+\\frac{t-t_i}{t_{i+1}-t_i}y(t_{i+1})
```
    
Constitute the original function with the interpolation and implement the summation.
"""
struct RLPiecewise <: RLInt end


"""
@book{oldham_spanier_1984,
title={The fractional calculus: Theory and applications of differentiation and integration to arbitrary order},
author={Oldham, Keith B. and Spanier, Jerome},
year={1984}}
"""

"""
# Riemann Liouville sense fractional integral approximation.

    fracint(f, α, end_point, h, RLIntApprox())

### Example

```julia-repl
julia> fracint(x->x^5, 0.5, 2.5, 0.0001, RLIntApprox())
64.36229646291393
```

Using the **Staircase approximation** to approximate the original function and implement the summation.
"""
struct RLIntApprox <: RLInt end


"""
# Riemann Liouville sense fractional integral using **Linear interpolation**.

    fracint(f, α, end_point, h, RLLinearInterp())

### Example

```julia-repl
julia> fracint(x->x^5, 0.5, 2.5, 0.0001, RLLinearInterp())
64.36234206434942
```

**RLLinearInterp** is more complex compared with *RLIntApprox* but more precise.

Deploying the [Linear Interpolation](https://en.wikipedia.org/wiki/Linear_interpolation) between ``f_{j+1}`` and ``f_j``, **RLLinearInterp** method is more precise than **RLIntApprox**.
"""
struct RLLinearInterp <: RLInt end


"""
@article{2009,
title={Matrix approach to discrete fractional calculus II: Partial fractional differential equations},
DOI={10.1016/j.jcp.2009.01.014},
author={Podlubny, Igor and Chechkin, Aleksei and Skovranek, Tomas and Chen, YangQuan and Vinagre Jara, Blas M.},
}

"""

"""
# Riemann Liouville sense integral using Triangular Strip Matrix to discrete.

    fracint(f, α, end_point, h, RLIntMatrix())

Using Triangular Strip Matrix to approximate fractional integral.

### Example

```julia-repl
julia> fracint(x->x^5, 0.5, 2.5, 0.0001, RLIntMatrix())
```

!!! info
    Triangular Strip Matrix method returns the derivative in the interval ``[0, T]`` in ```Vector```

!!! tip
    With the advancing Triangular Strip Matrix method, you can not only compute fractional integrals, integer order, higher order integral is also supported!!
Try to set α as an integer, arbitrary integer of course! I promise you would enjoy it😏

Using [Triangular Strip Matrix](https://en.wikipedia.org/wiki/Triangle_strip) to discrete the integral.
"""
struct RLIntMatrix <: RLInt end


"""
# Riemann Liouville sense fractional integral using fractional Simpson's formula to approximate.

### Example

```julia-repl
julia> fracint(x->x, 0.5, 1, 0.000001, RLIntSimpson())
0.7516516520541335
```

```tex
@inproceedings{Li2015NumericalMF,
  title={Numerical Methods for Fractional Calculus},
  author={Changpin Li and Fanhai Zeng},
  year={2015}
}
```

Using fractional Simpson's formula to discrete fractional integral.
"""
struct RLIntSimpson <: RLInt end


"""
# Riemann Liouville sense fractional integral using fractional trapezoidal formula to approximate.

### Example

```julia-repl
julia> fracint(x->x, 0.5, 1, 0.001, RLIntTrapezoidal())
0.7522527780636226
```

《Numerical methods for fractional calculus》.
Using Trapezoidal method to approximate fractional integral
"""
struct RLIntTrapezoidal <: RLInt end


"""
# Riemann Liouville sense fractional integral using fractional rectangular formula to approximate.

### Example

```julia-repl
julia> fracint(x->x, 0.5, 1, 0.001, RLIntRectangular())
0.7516812175993778
```
"""
struct RLIntRectangular <: RLInt end

"""
# Riemann Liouville sense fractional integral using cubic spline interpolation to approximate.

### Example

```julia-repl
julia> fracint(x->x, 0.5, 1, 0.0001, RLIntCubicSplineInterp())
0.7529119186452693
```

Error estimate of this method is ``\\mathcal{O(h^4)}``, it is determined by the error of the cubic spline interpolation.

!!! warning "Set h as 0.001 or bigger"
    For some reason, in the **RLIntCubicSplineInterp** method, set **h** as 0.001 would get better result.

> Numerical Methods for Fractional Calculus Page 34
"""
struct RLIntCubicSplineInterp <: RLInt end

####################################
###     Type definition done     ###
####################################



"""
Check if the format of nargins is correct
"""
function checks(f, α, start_point, end_point)
    α % 1 != 0 ? nothing : error("Decimal number only!")
    if isa(end_point, Number)
        start_point < end_point ? nothing : DomainError("Domain error! Start point must smaller than end point")
    elseif isa(end_point, AbstractArray)
        start_point < minimum(end_point) ? nothing : DomainError("Vector domain error! Start point must smaller than end point")
    else
        ErrorException("Please input correct point you wan tto compute")
    end

    #The fractional integral of number is relating with the end_point.
    if isa(f, Number)
        f == 0 ? 0 : 2*f*sqrt(end_point/pi)
    end
end

function fracint(f::FunctionAndNumber, α, start_point, end_point::Real, h::Real, ::RLDirect)
    #checks(f, α, start_point, end_point)
    typeof(f) <: Number ? (end_point == 0 ? (return 0) : (return 2*f*sqrt(end_point/pi))) : nothing
    end_point == 0 ? (return 0) : nothing

    temp1 = f(start_point) .* (end_point-start_point) .^α
    #Use Complex differentiation to obtain the differentiation
    g(τ) = imag(f(τ .+ 1*im*h) ./ h) .* (end_point-τ) .^α
    temp2 = quadgk(g, start_point, end_point)
    result = (temp1 .+ temp2) ./gamma(α+1)
    return result
end

function fracint(f::FunctionAndNumber, α::Real, end_point, h::Real, ::RLPiecewise)::Real
    #checks(f, α, 0, end_point)
    typeof(f) <: Number ? (end_point == 0 ? (return 0) : (return 2*f*sqrt(end_point/pi))) : nothing
    end_point == 0 ? (return 0) : nothing
    #Initialize
    n = round(Int, end_point/h)
    result = zero(Real)

    @fastmath @inbounds @simd for i ∈ 0:n
        result += W(i, n, α)*f(i*h)
    end

    result1 = result*h^α/gamma(α+2)
    return result1
end
function W(i, n, α)
    if i == 0
        return n^α*(α+1-n) + (n-1)^(α+1)
    elseif i == n 
        return 1
    else
        return (n-i-1)^(α+1) + (n-i+1)^(α+1) - 2*(n-i)^(α+1)
    end
end


function fracint(f::Union{Number, Function}, α::Real, end_point::AbstractArray, h::Real, ::RLPiecewise)::Vector
    result = map(x->fracint(f, α, x, h, RLPiecewise()), end_point)
    return result
end


function fracint(f::FunctionAndNumber, α::Real, end_point, h::Real, ::RLIntApprox)::Real
    #checks(f, α, 0, end_point)
    typeof(f) <: Number ? (end_point == 0 ? (return 0) : (return 2*f*sqrt(end_point/pi))) : nothing
    end_point == 0 ? (return 0) : nothing
    α = -α
    n = round(Int, end_point/h)
    result = zero(Real)

    @fastmath @inbounds @simd for i ∈ 0:n-1
        result += (f(end_point-i*h) + f(end_point-(i+1)*h))*((i+1)^(-α) - i^(-α))
    end

    result1 = result/2*end_point^(-α)*n^α/gamma(1-α)
    return result1
end

function fracint(f::Union{Number, Function}, α::Real, end_point::AbstractArray, h::Real, ::RLIntApprox)::Vector
    result = map(x->fracint(f, α, x, h, RLIntApprox()), end_point)
    return result
end




function fracint(f::FunctionAndNumber, α::Real, end_point::Number, h::Real, ::RLLinearInterp)::Real
    typeof(f) <: Number ? (end_point == 0 ? (return 0) : (return 2*f*sqrt(end_point/pi))) : nothing
    end_point == 0 ? (return 0) : nothing

    α = -α
    n = round(Int, end_point/h)
    result = zero(Real)
    
    @fastmath @inbounds @simd for i ∈ 0:n-1
        result += ((i+1)*f(end_point-i*h)-i*f(end_point-(i+1)*h))/(-α)*((i+1)^(-α) - i^(-α))+(f(end_point - (i+1)*h) - f(end_point-i*h))/(1-α)*((i+1)^(1-α) - i^(1-α))
    end

    result1 = result*end_point^(-α)*n^α/gamma(-α)
    return result1
end

function fracint(f::Union{Number, Function}, α::Real, end_point::AbstractArray, h::Real, ::RLLinearInterp)::Vector
    result = map(x->fracint(f, α, x, h, RLLinearInterp()), end_point)
    return result
end



function fracint(f, α::Number, end_point, h::Real, ::RLIntMatrix)
    N = round(Int, end_point/h+1)
    tspan = collect(0:h:end_point)
    return J(N, α, h)*f.(tspan)
end

function ω(n, p)
    omega = zeros(n+1)
    omega[1]=1
    @fastmath @inbounds @simd for i ∈ 1:n
        omega[i+1]=(1-(p+1)/i)*omega[i]
    end
    
    return omega
end

function J(N, p, h::Real)
    result = zeros(N, N)
    temp = ω(N, -p)

    @fastmath @inbounds @simd for i ∈ 1:N
        result[i, 1:i] = reverse(temp[1:i])
    end

    return h^p*result
end


#=
RLIntSimpson Algorithm
=#
function fracint(f, α, point, h::Real, ::RLIntSimpson)
    typeof(f) <: Number ? (point == 0 ? (return 0) : (return 2*f*sqrt(point/pi))) : nothing
    point == 0 ? (return 0) : nothing

    N = round(Int, point/h)
    temp1 = 0
    temp2 = 0

    @fastmath @inbounds @simd for i ∈ 0:N
        temp1 += cₖₙ(i, N, α)*f(i*h)
    end

    @fastmath @inbounds @simd for k ∈ 0:N-1
        temp2 += ĉₖₙ(k, N, α)*f((k+0.5)*h)
    end

    return h^α/gamma(α+3)*temp1 + 4*h^α/gamma(α+3)*temp2
end

function cₖₙ(k, n, α)
    if k == 0
        return 4*((n+1)^(2+α)-n^(2+α))-(α+2)*(3*(n+1)^(1+α)+n^(1+α))+(α+2)*(α+1)*(n+1)^α
    elseif 1 ≤ k ≤ n-1
        return -(α+2)*((n+1-k)^(1+α)+6*(n-k)^(1+α)+(n-k-1)^(1+α))+4*((n+1-k)^(2+α)-(n-1-k)^(2+α))
    elseif k == n
        return 2-α
    end
end
function ĉₖₙ(k, n, α)
    return ((α+2)*((n+1-k)^(1+α)+(n-k)^(1+α))-2*((n+1-k)^(2+α)-(n-k)^(2+α)))
end

function fracint(f::Union{Number, Function}, α::Real, end_point::AbstractArray, h::Real, ::RLIntSimpson)::Vector
    result = map(x->fracint(f, α, x, h, RLIntSimpson()), end_point)
    return result
end

#=
RLIntTrapezoidal Algorithm
=#
function fracint(f::FunctionAndNumber, α, point, h::Real, ::RLIntTrapezoidal)
    typeof(f) <: Number ? (point == 0 ? (return 0) : (return 2*f*sqrt(point/pi))) : nothing
    point == 0 ? (return 0) : nothing

    N = round(Int, point/h)
    result = zero(Real)

    @fastmath @inbounds @simd for i ∈ 0:N
        result += aₖₙ(i, N, α)*f(i*h)
    end

    return h^α/gamma(α+2)*result
end

function aₖₙ(k, n, α)
    if k == 0
        return (n-1)^(α+1) - (n-1-α)*n^α
    elseif 1 ≤ k ≤ n-1
        return (n-k+1)^(α+1)+(n-1-k)^(α+1)-2*(n-k)^(α+1)
    else
        return 1
    end
end


function fracint(f::Union{Number, Function}, α, point, h::Real, ::RLIntRectangular)
    typeof(f) <: Number ? (point == 0 ? (return 0) : (return 2*f*sqrt(point/pi))) : nothing
    point == 0 ? (return 0) : nothing

    N = round(Int, point/h)
    result = zero(Real)

    @fastmath @inbounds @simd for i ∈ 0:N-1
        result += bₖ(N-i-1, α)*f(i*h)
    end

    return h^α/gamma(α+1)*result
end

function bₖ(k, α)
    return (k+1)^α-k^α
end

function fracint(f::Union{Number, Function}, α::Real, end_point::AbstractArray, h::Real, ::RLIntRectangular)::Vector
    result = map(x->fracint(f, α, x, h, RLIntRectangular()), end_point)
    return result
end

#=
RLIntCubicSplineInterp Algorithm, when h is 0.01 behave best
=#
function fracint(f, α, point, h::Real, ::RLIntCubicSplineInterp)
    typeof(f) <: Number ? (point == 0 ? (return 0) : (return 2*f*sqrt(point/pi))) : nothing
    point == 0 ? (return 0) : nothing
    
    N = round(Int, point/h)
    result = zero(Real)

    @fastmath @inbounds @simd for j ∈ 0:N
        result += eⱼₙ(j, N, α)*f(j*h) + h*êⱼₙ(j, N, α)*first_order(f, j*h, h)
    end

    return h^α/gamma(α+4)*result
end

function eⱼₙ(j, n, α)
    if j == 0
        return -6*(n-1)^(2+α)*(1+2*n+α) + n^α*(12*n^3-6*(3+α)*n^2 + (1+α)*(2+α)*(3+α))
    elseif j == n
        return 6*(1+α)
    else
        return 6*(4*(n-j)^(3+α) + (n-j-1)^(2+α)*(2*j-2*n-1-α) + (1+n-j)^(2+α)*(2*j-2*n+1+α))
    end
end
function êⱼₙ(j, n, α)
    if j == 0
        return -2*(n-1)^(2+α)*(3*n+α) + n^(1+α)*(6*n^2-4*(3+α)*n + (2+α)*(3+α))
    elseif j == n
        return -2*α
    else
        return 2*(n-j-1)^(α+2)*(3*j-3*n-α)-2*(n-j+1)^(α+2)*(3*j-3*n+α) - (n-j)^(α+2)*(24+8α)
    end
end

function fracint(f::Union{Number, Function}, α::Real, end_point::AbstractArray, h::Real, ::RLIntCubicSplineInterp)::Vector
    result = map(x->fracint(f, α, x, h, RLIntCubicSplineInterp()), end_point)
    return result
end