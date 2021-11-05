"""
The base type of a fractional integral sense, in FractionalCalculus.jl, all of the fractional integral senses belong to ```FracIntAlg```
"""
abstract type FracIntAlg end

"""
Riemann-Liouville sense fractional integral algorithms

Note this two algorithms belong to direct compute, precise are ensured, but maybe cause more memory allocation and take more compilation time.
"""
abstract type RLInt <: FracIntAlg end

"""
Using the direct mathematic expression:

```math
(J^αf)(t)=\\frac{1}{\\Gamma(α)} \\int_0^t(t-τ)^{α-1}f(τ)dτ
```

By using [QuadGK](https://github.com/JuliaMath/QuadGK.jl) calculate the integration and obtain the value.
"""
struct RL_Direct <: RLInt end

"""
With first derivative known, we can use the Riemann Liouville sense to obtain the fractional integral more effcient.
"""
struct RL_Direct_First_Diff_Known <: RLInt end

"""
@article{LI20113352,
title = {Numerical approaches to fractional calculus and fractional ordinary differential equation},
author = {Changpin Li and An Chen and Junjie Ye},
}
"""

"""
Using piecewise linear interpolation:

```math
    y_n(t)=\\frac{t-t_{i+1}}{t_i-t_{i+1}}y(t_i)+\\frac{t-t_i}{t_{i+1}-t_i}y(t_{i+1})
```

Constitute the original function with the interpolation and implement the summation.
"""
struct RL_Piecewise <: RLInt end


"""
@book{oldham_spanier_1984,
title={The fractional calculus: Theory and applications of differentiation and integration to arbitrary order},
author={Oldham, Keith B. and Spanier, Jerome},
year={1984}}
"""

"""
Using the **Staircase approximation** to approximate the original function and implement the summation.
"""
struct RLInt_Approx <: RLInt end

"""
Deploying the [Linear Interpolation](https://en.wikipedia.org/wiki/Linear_interpolation) between ``f_{j+1}`` and ``f_j``, **RL_LinearInterp** method is more precise than **RLInt_Approx**.
"""
struct RL_LinearInterp <: RLInt end


"""
@article{2009,
title={Matrix approach to discrete fractional calculus II: Partial fractional differential equations},
DOI={10.1016/j.jcp.2009.01.014},
author={Podlubny, Igor and Chechkin, Aleksei and Skovranek, Tomas and Chen, YangQuan and Vinagre Jara, Blas M.},
}

"""

"""
Using [Triangular Strip Matrix](https://en.wikipedia.org/wiki/Triangle_strip) to discrete the integral.
"""
struct RLInt_Matrix <: RLInt end


####################################
###     Type defination done     ###
####################################



"""
Check if the format of nargins is correct
"""
function checks(α, start_point, end_point)
    α % 1 != 0 ? nothing : error("Decimal number only!")
    if typeof(end_point) <: Number
        start_point < end_point ? nothing : DomainError("Domain error! Start point must smaller than end point")
    elseif typeof(end_point) <: AbstractArray
        start_point < minimum(end_point) ? nothing : DomainError("Vector domain error! Start point must smaller than end point")
    else
        ErrorException("Please input correct point you wan tto compute")
    end   
end


"""
# Riemann Liouville sense fractional integral

    fracint(f::Function, α, start_point, end_point, h, RL_Direct())

### Example: 

```julia-repl
julia> fracint(x->x^5, 0.5, 0, 2.5, 1e-8, RL_Direct())
```

Riemann_Liouville fractional integral using complex step differentiation.
Returns a tuple (1.1639316474512205, 1.0183453796725215e-8), which contains the value of this derivative is 1.1639316474512205, and the error estimate is 1.0183453796725215e-8
"""
function fracint(f::Union{Function, Number}, α, start_point, end_point::Real, h, ::RL_Direct)
    checks(α, start_point, end_point)
    
    #The fractional integral of number is relating with the end_point.
    if typeof(f) <: Number
        if f == 0
            return 0
        else
            return 2*f*sqrt(end_point/pi)
        end
    end

    temp1 = f(start_point) .* (end_point-start_point) .^α
    #Use Complex differentiation to obtain the differentiation
    g(τ) = imag(f(τ .+ 1*im*h) ./ h) .* (end_point-τ) .^α
    temp2 = quadgk(g, start_point, end_point)
    result = (temp1 .+ temp2) ./gamma(α+1)
    return result
end

function fracint(f::Union{Function, Number}, α::Float64, start_point, end_point::AbstractArray, h, ::RL_Direct)::Vector
    ResultArray = Float64[]

    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracint(f, α, start_point, value, h, RL_Direct())[1])
    end
    return ResultArray
end


"""
# Riemann Liouville sense fractional integral with first diff known.

    fracint(f, fd, α, start_point, end_point, RL_Direct_First_Diff_Known())

### Example

```julia-repl
julia> fracint(x->x^5, x->5*x^4, 0.5, 0, 2.5, RL_Direct_First_Diff_Known())
```

With first order derivative known, we can directly use it in the computation of α order fraction integral
"""
function fracint(f::Function, fd::Function, α, start_point, end_point, ::RL_Direct_First_Diff_Known)
    checks(α, start_point, end_point)
    
    #The fractional integral of number is relating with the end_point.
    if typeof(f) <: Number
        if f == 0
            return 0
        else
            return 2*f*sqrt(end_point/pi)
        end
    end
    
    temp1 = f(start_point) .* (end_point-start_point) .^α
    g(τ) = fd(τ) .* (end_point-τ) .^α
    temp2 = quadgk(g, start_point, end_point)
    result = (temp1 .+ temp2) ./ gamma(α+1)
    return result
end

function fracint(f::Union{Function, Number}, α::Float64, start_point, end_point::AbstractArray, ::RL_Direct_First_Diff_Known)::Vector
    ResultArray = Float64[]

    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracint(f, α, start_point, value, RL_Direct_First_Diff_Known())[1])
    end
    return ResultArray
end


"""
# Riemann Liouville sense fractional integral using piecewise interpolation.

    fracint(f, α, end_point, h, RL_Piecewise())

### Example

```julia-repl
julia> fracint(x->x^5, 0.5, 2.5, 0.0001, RL_Piecewise())
```

By deploying Piecewise interpolation to approximate the original function, with small step size, this method is fast and take little memory allocation.
"""
function fracint(f::Union{Function, Number}, α::Float64, end_point, h, ::RL_Piecewise)::Float64
    end_point > 0 ? nothing : error("Please compute the integral of a positive value")
    
    #The fractional integral of number is relating with the end_point.
    if typeof(f) <: Number
        if f == 0
            return 0
        else
            return 2*f*sqrt(end_point/pi)
        end
    end

    n=end_point/h
    result=0

    for i in range(0, n, step=1)
        result = result + W(i, n, α)*f(i*h)
    end

    result1 = result*h^α/gamma(α+2)
    return result1
end
function W(i, n, α)
    if i==0
        return n^α*(α+1-n)+(n-1)^(α+1)
    elseif i==n 
        return 1
    else
        return (n-i-1)^(α+1)+(n-i+1)^(α+1)-2*(n-i)^(α+1)
    end
end


function fracint(f::Union{Function, Number}, α::Float64, end_point::AbstractArray, h, ::RL_Piecewise)::Vector
    ResultArray = Float64[]

    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracint(f, α, value, h, RL_Piecewise())[1])
    end
    return ResultArray
end


"""
# Riemann Liouville sense fractional integral approximation.

    fracint(f, α, end_point, h, RLInt_Approx())

### Example

```julia-repl
julia> fracint(x->x^5, 0.5, 2.5, 0.0001, RLInt_Approx())
```
"""
function fracint(f::Union{Function, Number}, α::Float64, end_point, h, ::RLInt_Approx)::Float64
        
    #The fractional integral of number is relating with the end_point.
    if typeof(f) <: Number
        if f == 0
            return 0
        else
            return 2*f*sqrt(end_point/pi)
        end
    end

    if end_point == 0
        return 0
    end

    α = -α
    n = end_point/h
    result=0

    for i in range(0, n-1, step=1)
        result += (f(end_point-i*end_point/n)+f(end_point-(i+1)*end_point/n))/2*((i+1)^(-α)-i^(-α))
    end

    result1=result*end_point^(-α)*n^α/gamma(1-α)
    return result1
end

function fracint(f::Union{Function, Number}, α::Float64, end_point::AbstractArray, h, ::RLInt_Approx)::Vector
    ResultArray = Float64[]

    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracint(f, α, value, h, RLInt_Approx())[1])
    end
    return ResultArray
end


"""
# Riemann Liouville sense fractional integral using Linear interpolation.

    fracint(f, α, end_point, h, RL_LinearInterp())

### Example

```julia-repl
julia> fracint(x->x^5, 0.5, 2.5, 0.0001, RL_LinearInterp())
```

**RL_LinearInterp** is more complex but more precise.
"""
function fracint(f::Union{Function, Number}, α::Float64, end_point, h, ::RL_LinearInterp)::Float64
        
    #The fractional integral of number is relating with the end_point.
    if typeof(f) <: Number
        if f == 0
            return 0
        else
            return 2*f*sqrt(end_point/pi)
        end
    end

    if end_point == 0
        return 0
    end

    #Support both Matrix and Vector end points
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            append!(ResultArray, fracint(f, α, value, h, RL_LinearInterp()))
        end
        return ResultArray
    end

    α = -α
    n = end_point/h
    result = 0
    
    for i in range(0, n-1, step=1)
        result+=((i+1)*f(end_point-i*end_point/n)-i*f(end_point-(i+1)*end_point/n))/(-α)*((i+1)^(-α)-i^(-α))+(f(end_point-(i+1)*end_point/n)-f(end_point-i*end_point/n))/(1-α)*((i+1)^(1-α)-i^(1-α))
    end

    result1=result*end_point^(-α)*n^α/gamma(-α)
    return result1
end

function fracint(f::Union{Function, Number}, α::Float64, end_point::AbstractArray, h, ::RL_LinearInterp)::Vector
    ResultArray = Float64[]

    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracint(f, α, value, h, RL_LinearInterp())[1])
    end
    return ResultArray
end

"""
# Riemann Liouville sense integral using Triangular Strip Matrix to discrete and compute.

    fracint(f, α, end_point, h, RLInt_Matrix())

Using Triangular Strip Matrix to approximate fractional integral.

### Example

```julia-repl
julia> fracint(x->x^5, 0.5, 2.5, 0.0001, RLInt_Matrix())
```

!!! info
    Triangular Strip Matrix method returns the derivative in the interval ``[0, T]`` in ```Vector```

!!! tip
    With the advancing Triangular Strip Matrix method, you can not only compute fractional integrals, integer order, higher order integral is also supported!!
Try to set α as an integer, arbitrary integer of course! I promise you would enjoy it😏

"""
function fracint(f, α, end_point, h, ::RLInt_Matrix)
    N=Int64(end_point/h+1)
    tspan=collect(0:h:end_point)
    return J(N, α, h)*f.(tspan)
end

function omega(n, p)
    omega = zeros(n+1)

    omega[1]=1
    for i in range(1, n, step=1)
        omega[i+1]=(1-(p+1)/i)*omega[i]
    end
    
    return omega

end
function J(N, p, h)
    result=zeros(N, N)
    temp=omega(N, -p)

    for i in range(1, N, step=1)
        result[i, 1:i]=reverse(temp[1:i])
    end

    return h^p*result
end




## Macros for convenient computing.
"""
    @fracint(f, α, point)

Return the α-order integral of **f** at specific point.

```julia-repl
julia> @fracint(x->x, 0.5, 1)
0.7522525439593486
```
"""
macro fracint(f::Union{Function, Number}, α::Float64, point)
    return :(fracint($f, $α, $point, 0.0001, RLInt_Approx()))
end

"""
    @semifracint(f, point)

Return the semi-integral of **f** at specific point.

```julia-repl
julia> @semifracint(x->x, 1)
0.7522525439593486
```
"""
macro semifracint(f::Union{Function, Number}, point)
    return :(fracint($f, 0.5, $point, 0.0001, RLInt_Approx()))
end