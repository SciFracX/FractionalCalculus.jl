import FractionalCalculus.FracDiffAlg

"""
Riemann-Liouville sense fractional derivative algorithms, please refer to [Riemann-Liouville derivative](https://en.wikipedia.org/wiki/Fractional_calculus#Riemann%E2%80%93Liouville_fractional_derivative)
"""
abstract type RLDiff <: FracDiffAlg end

#RLDiff_LinearInterp maybe??
"""
Using Linear interpolation to approximate fractional derivative in Riemann Liouville  fractional derivative sense.
"""
struct RLDiff_Approx <: RLDiff end

"""
@article{2009,
title={Matrix approach to discrete fractional calculus II: Partial fractional differential equations},
DOI={10.1016/j.jcp.2009.01.014},
author={Podlubny, Igor and Chechkin, Aleksei and Skovranek, Tomas and Chen, YangQuan and Vinagre Jara, Blas M.},
}

"""

"""
Using [Triangular Strip Matrix](https://en.wikipedia.org/wiki/Triangle_strip) to discrete the derivative.
"""
struct RLDiff_Matrix <: RLDiff end


"""
Using linear spline interpolation method to approximate the Riemann Liouville fractional derivative.
"""
struct RL_Linear_Spline_Interp <: RLDiff end


################################################################
###                    Type defination done                  ###
################################################################


"""
# Riemann Liouville sense derivative approximation

    fracdiff(f, α, end_point, h, RLDiff_Approx())

Using approximation to obtain fractional derivative value.

### Example

```julia-repl
julia> fracdiff(x->x^5, 0.5, 2.5, 0.0001, RLDiff_Approx())
```

!!! warning
    The RLDiff_Approx algorithm only support for 0 < α < 1.

"""
function fracdiff(f::Union{Number, Function}, α, end_point, h, ::RLDiff_Approx)
    #checks(f, α, 0, end_point)

    summation = 0
    n = Int64(floor(end_point/h))

    @fastmath @inbounds @simd for i ∈ 0:n-1
        summation += (f(end_point-i*h) - f(end_point-(i+1)*h))*((i+1)^(1-α) - i^(1-α))
    end

    result = ((1-α)*f(0)/n^α+summation)*end_point^(-α)*n^α/gamma(2-α)
    return result
end

function fracdiff(f::Union{Number, Function}, α, end_point::AbstractArray, h, ::RLDiff_Approx)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, value, h, RLDiff_Approx()))
    end
    return ResultArray
end


"""
# Riemann Liouville sense derivative using Triangular Strip Matrix to discrete and compute.

    fracdiff(f, α, end_point, h, RLDiff_Matrix())

Using Triangular Strip Matrix to approximate fractional derivative.

### Example

```julia-repl
julia> fracdiff(x->x^5, 0.5, 2.5, 0.0001, RLInt_Matrix())
```

!!! info
    Triangular Strip Matrix method returns the derivative in the interval ``[0, T]`` in ```Vector```

!!! tip
    With the advancing Triangular Strip Matrix method, you can not only compute fractional derivatives, integer order, higher order derivative is also supported!!
Try to set α as an integer, arbitrary integer of course! I promise you would enjoy it😏
"""
function fracdiff(f, α, end_point, h::Float64, ::RLDiff_Matrix)
    N = Int64(end_point/h+1)
    @views tspan = collect(0:h:end_point)
    return B(N, α, h)*f.(tspan)
end

#Compute the eliminator matrix Sₖ by omiting n-th row
function eliminator(n, row)
    temp = zeros(n, n) + I
    return @views temp[Not(row), :]
end

#Generate elements in Matrix.
function omega(n, p)
    omega = zeros(n+1)

    omega[1]=1
    @fastmath @inbounds @simd for i ∈ 1:n
        omega[i+1]=(1-(p+1)/i)*omega[i]
    end
    
    return omega
end
function B(N, p, h::Float64)
    result = zeros(N, N)
    temp = omega(N, p)

    @inbounds @simd for i ∈ 1:N
        @views result[i, 1:i]=reverse(temp[1:i])
    end

    return h^(-p)*result
end
# Multiple dispatch for not assigning step size *h*.
function B(N, p)
    result=zeros(N, N)
    temp=omega(N, p)

    @inbounds @simd for i ∈ 1:N
        @views result[i, 1:i]=reverse(temp[1:i])
    end

    return result
end



"""
Numerical methods for fractional calculus by Li, Changpin
Page 57

Linear Spline Interpolation
"""
function fracdiff(f, α, x, h, ::RL_Linear_Spline_Interp)
    N = Int64(floor(x/h))

    result = 0

    for k = 0:(N+1)
        result += z̄ₘₖ(N, k, α)*f(k*h)
    end

    return 1/(gamma(4-α)h^α)*result

end

function z̄ₘₖ(m, k, α)
    if k ≤ m-1
        return c̄ⱼₖ(m-1, k, α)-2*c̄ⱼₖ(m, k, α)+c̄ⱼₖ(m+1, k, α)
    elseif k == m
        return -2*c̄ⱼₖ(m, k, α)+c̄ⱼₖ(m+1, k, α)
    elseif k == m+1
        return c̄ⱼₖ(m+1, k, α)
    elseif k > m+1
        return 0
    end
end

function c̄ⱼₖ(j, k, α)
    if k==0
        return (j-1)^(3-α)-j^(2-α)*(j-3+α)
    elseif 1 ≤ k ≤ j-1
        return (j-k+1)^(3-α)-2*(j-k)^(3-α)+(j-k-1)^(3-α)
    elseif k == j
        return 1
    end
end