"""
The base type of a fractional integral sense, in FractionalCalculus.jl, all of the fractional integral senses belong to ```FracIntAlg```
"""
abstract type FracIntAlg end

"""
Riemann-Liouville sense fractional integral algorithms

Note this two algorithms belong to direct compute, precise are ensured, but maybe cause more memory allocation and take more compilation time.
"""
abstract type RL <: FracIntAlg end
struct RL_Direct <: RL end
struct RL_Direct_First_Diff_Known <: RL end

"""
@article{LI20113352,
title = {Numerical approaches to fractional calculus and fractional ordinary differential equation},
author = {Changpin Li and An Chen and Junjie Ye},
}

Using piecewise linear interpolation function to approximate input function and implement summation.
"""
struct RL_Piecewise <: RL end



"""
Check if the format of margins is correct
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
**fracint** as the general API used to compute fractional integral.
By entering the function, the order, the start point and end point and step size, and choose the algorithms you want to use
    fracint(f::Function)

# Example: 

```julia
fracint(x->x^5, 0.5, 0, 2.5, 1e-8, RL_Direct())
```

Riemann_Liouville fractional integral using complex step differentiation
Returns a tuple (1.1639316474512205, 1.0183453796725215e-8), which means the value of this derivative is 1.1639316474512205, and the error estimate is 1.0183453796725215e-8
"""
function fracint(f::Union{Function, Number}, α, start_point, end_point, step_size, ::RL_Direct)
    checks(α, start_point, end_point)
    
    #Support vectorized end point
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            #println(fracdiff(f, α, start_point, value, step_size, RL()))
            append!(ResultArray, fracint(f, α, start_point, value, step_size, RL_Direct())[1])
        end
        return ResultArray
    end

    temp1 = f(start_point) .* (end_point-start_point) .^α
    #Use Complex differentiation to obtain the differentiation
    g(τ) = imag(f(τ .+ 1*im*step_size) ./ step_size) .* (end_point-τ) .^α
    temp2 = quadgk(g, start_point, end_point)
    result = (temp1 .+ temp2) ./gamma(α+1)
    return result
end


"""
With first order derivative known, we can directly use it in the computation of α order fraction integral
"""
function fracint(f::Function, fd::Function, α, start_point, end_point, ::RL_Direct_First_Diff_Known)
    checks(α, start_point, end_point)
    
    #Support vectorized end point
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            #println(fracdiff(f, α, start_point, value, step_size, RL()))
            append!(ResultArray, fracint(f, α, start_point, value, step_size, RL_Direct_First_Diff_Known())[1])
        end
        return ResultArray
    end
    
    temp1 = f(start_point) .* (end_point-start_point) .^α
    g(τ) = fd(τ) .* (end_point-τ) .^α
    temp2 = quadgk(g, start_point, end_point)
    result = (temp1 .+ temp2) ./ gamma(α+1)
    return result
end

"""
By deploying Piecewise interpolation to approximate the original function, with small step_size, this method is fast and take little memory allocation.
"""
function fracint(f::Function, α, end_point, step_size, ::RL_Piecewise)
    end_point > 0 ? nothing : error("Please compute the integral of a positive value")
    
    #Support vectorized end point
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            append!(ResultArray, fracint(f, α, value, step_size, RL_Piecewise()))
        end
        return ResultArray
    end

    n=end_point/step_size
    result=0

    for i in range(0, n, step=1)
        result = result + W(i, n, α)*f(i*step_size)
    end

    result1 = result*step_size^α/gamma(α+2)
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