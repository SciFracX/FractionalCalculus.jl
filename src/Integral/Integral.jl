"""
The base type of a fractional integral sense, in FractionalCalculus.jl, all of the fractional integral senses belong to ```FracIntAlg```
"""
abstract type FracIntAlg end

"""
Riemann-Liouville sense fractional integral algorithms
"""
struct RL <: FracIntAlg end
struct RL_First_Diff_Known <: FracIntAlg end



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
fracdiff(x->x^5, 0.5, 0, 2.5, 1e-8)
```

Riemann_Liouville fractional integral using complex step differentiation
Returns a tuple (1.1639316474512205, 1.0183453796725215e-8), which means the value of this derivative is 1.1639316474512205, and the error estimate is 1.0183453796725215e-8
"""
function fracint(f, α, start_point, end_point, step_size, ::RL)
    checks(α, start_point, end_point)
    
    #Support vectorized end point
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            #println(fracdiff(f, α, start_point, value, step_size, RL()))
            append!(ResultArray, fracint(f, α, start_point, value, step_size, RL())[1])
        end
        return ResultArray
    end

    temp1 = f(start_point) .* (end_point-start_point) .^α
    g(τ) = imag(f(τ .+ 1*im*step_size) ./ step_size) .* (end_point-τ) .^α
    temp2 = quadgk(g, start_point, end_point)
    result = (temp1 .+ temp2) ./gamma(α+1)
    return result
end


"""
With first order derivative known, we can directly use it in the computation of α order fraction integral
"""
function fracint(f::Function, fd::Function, α, start_point, end_point, ::RL_First_Diff_Known)
    checks(α, start_point, end_point)
    
    #Support vectorized end point
    if typeof(end_point) <: AbstractArray
        ResultArray = Float64[]
        for (_, value) in enumerate(end_point)
            #println(fracdiff(f, α, start_point, value, step_size, RL()))
            append!(ResultArray, fracint(f, α, start_point, value, step_size, RL_First_Diff_Known())[1])
        end
        return ResultArray
    end
    
    temp1 = f(start_point) .* (end_point-start_point) .^α
    g(τ) = fd(τ) .* (end_point-τ) .^α
    temp2 = quadgk(g, start_point, end_point)
    result = (temp1 .+ temp2) ./ gamma(α+1)
    return result
end