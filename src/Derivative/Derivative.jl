"""
Base type of fractional differentiation senses, , in FractionalCalculus.jl, all of the fractional derivative senses belong to ```FracDiffAlg```
"""
abstract type FracDiffAlg end
abstract type Caputo <: FracDiffAlg end
abstract type GL <: FracDiffAlg end

"""
Note these four algorithms belong to direct compute, precise are ensured, but maybe cause more memory allocation and take more compilation time.
"""
struct Caputo_Direct <: Caputo end
struct Caputo_Direct_First_Diff_known <: Caputo end
struct Caputo_Direct_First_Second_Diff_Known <: Caputo end

struct GL_Direct <: GL end


"""
@article{LI20113352,
title = {Numerical approaches to fractional calculus and fractional ordinary differential equation},
author = {Changpin Li and An Chen and Junjie Ye},
}

Using piecewise linear interpolation function to approximate input function and combining Caputo derivative then implement summation.
"""
struct Caputo_Piecewise <: Caputo end


"""
Check if the format of margins is correct
"""
function checks(α, start_point, end_point)
    α % 1 != 0 ? nothing : error("α must be a decimal number")
    if typeof(end_point) <: Number
        start_point < end_point ? nothing : DomainError("Domain error! Start point must smaller than end point")
    elseif typeof(end_point) <: AbstractArray
        start_point < minimum(end_point) ? nothing : DomainError("Vector domain error! Start point must smaller than end point")
    else
        ErrorException("Please input correct point you wan tto compute")
    end
end


"""
    fracdiff(f::Function, α, start_point, end_point, ::FracDiffAlg)

Compute the α-order fractional derivative of f from start point to end point with specific algorithm.
"""
function fracdiff(f, α, start_point, end_point, ::FracDiffAlg)
end


"""
    fracdiff(f::Function, α, start_point, end_point, step_size ::Caputo)

Using Caputo definition to compute fractional derivative with **complex step differentiation**.

# Example: 

```julia
fracdiff(x->x^5, 0.5, 0, 2.5, 0.0001, Caputo_Direct())
```

Returns a tuple (**result**, **error**), which means the value of this derivative is 141.59714979764541, and the error estimate is 1.1532243848672914e-6.

When the input end points is an array, **fracdiff** will compute 

https://en.wikipedia.org/wiki/Fractional_calculus#Caputo_fractional_derivative

"""
function fracdiff(f, α, start_point, end_point, step_size, ::Caputo_Direct)
    checks(α, start_point, end_point)

    #The fractional derivative of a constant is zero
    if typeof(f) <: Number
        return 0
    end
    
    #Support both Matrix and Vector end points
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            #println(fracdiff(f, α, start_point, value, step_size, Caputo()))
            append!(ResultArray, fracdiff(f, α, start_point, value, step_size, Caputo_Direct())[1])
        end
        return ResultArray
    end

    if end_point == 0
        return 0
    end

    #Using complex step differentiation to calculate the first order differentiation
    g(τ) = imag(f(τ+1*im*step_size) ./ step_size) ./ ((end_point-τ) .^α)
    result = quadgk(g, start_point, end_point) ./gamma(1-α)
    return result
end


"""
    fracdiff(f, α, start_point, end_point, ::GL)

Using Grunwald-Letnikov definition to compute fractional derivative.

# Example:

```julia
fracdiff(x->x^5, 0, 0.5, GL_Direct())
```

> Please note Grunwald-Letnikov sense fracdiff only support 0 < α < 1.

Please refer to https://en.wikipedia.org/wiki/Gr%C3%BCnwald%E2%80%93Letnikov_derivative.
"""
function fracdiff(f::Union{Function, Number}, α, start_point, end_point, ::GL_Direct)
    checks(α, start_point, end_point)

    #The fractional derivative of a constant is zero
    if typeof(f) <: Number
        return 0
    end
    
    #Support both Matrix and Vector end points
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            append!(ResultArray, fracdiff(f, α, start_point, value, GL_Direct())[1])
        end
        return ResultArray
    end

    checks(α, start_point, end_point)

    if end_point == 0
        return 0
    end

    g(τ) = (f(end_point)-f(τ)) ./ (end_point - τ) .^ (1+α)
    result = f(end_point)/(gamma(1-α) .* (end_point - start_point) .^ α) .+ quadgk(g, start_point, end_point) .* α ./ gamma(1-α)
    return result
end


"""
If the first order derivative of a function is already known,

We can use this method to compute the fractional order derivative more precisely.

"""
function fracdiff(fd::Function, α, start_point, end_point, ::Caputo_Direct_First_Diff_known)
    checks(α, start_point, end_point)

    #The fractional derivative of a constant is zero
    if typeof(f) <: Number
        return 0
    end
    
    #Support both Matrix and Vector end points
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            #println(fracdiff(f, α, start_point, value, step_size, Caputo()))
            append!(ResultArray, fracdiff(f, α, start_point, value, Caputo_Direct_First_Diff_known())[1])
        end
        return ResultArray
    end

    g(τ) = (end_point-τ) .^ (-α) * fd(τ)
    result = quadgk(g, start_point, end_point) ./gamma(1-α)    
    return result
end


"""
If the first and second order derivative of a function is already known
"""
function fracdiff(fd1::Function, fd2, α, start_point, end_point, ::Caputo_Direct_First_Second_Diff_Known)
    checks(α, start_point, end_point)
    
    #Support both Matrix and Vector end points
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            append!(ResultArray, fracdiff(f, α, start_point, value, step_size, Caputo_Direct_First_Second_Diff_Known())[1])
        end
        return ResultArray
    end

    temp1 = fd1(start_point) .* (end_point-start_point) .^ (1-α) ./ gamma(2-α)
    g(τ) = fd2(τ) .* ((end_point-τ) .^ (1-α))
    temp2 = quadgk(g, start_point, end_point) ./ gamma(2-α)
    result = temp1 .+ temp2
    return result
end

"""
Using the piecewise algorithm to obtain the fractional derivative at a specific point.
"""
function fracdiff(f, α, end_point, step_size, ::Caputo_Piecewise)

    if typeof(f) <: Number
        return 0
    end

    if end_point == 0
        return 0
    end

    #Support both Matrix and Vector end points
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            append!(ResultArray, fracdiff(f, α, value, step_size, Caputo_Piecewise()))
        end
        return ResultArray
    end

    m = floor(α)+1

    result = 0
    n = end_point/step_size

    for i in range(0, n, step=1)
        result +=W₅(i, n, m, α)*first_order(f, i*step_size, step_size)
    end

    result1=result*step_size^(m-α)/gamma(m-α+2)
    return result1
end
function W₅(i, n, m, α)
    if i==0
        return n^(m-α)*(m-α+1-n)+(n-1)^(m-α+1)
    elseif i==n
        return 1
    else
        (n-i-1)^(m-α+1)+(n-i+1)^(m-α+1)-2*(n-i)^(m-α+1)
    end
end
function first_order(f, point, h)
    return (f(point+h)-f(point-h))/(2*h)
end