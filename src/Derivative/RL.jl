import FractionalCalculus.FracDiffAlg

"""
Riemann-Liouville sense fractional derivative algorithms, please refer to [Riemann-Liouville derivative](https://en.wikipedia.org/wiki/Fractional_calculus#Riemann%E2%80%93Liouville_fractional_derivative)
"""
abstract type RLDiff <: FracDiffAlg end

#RLDiff_LinearInterp maybe??
"""
# Riemann Liouville sense derivative approximation

    fracdiff(f, α, end_point, h, RLDiff_Approx())

Using Linear interpolation to approximate fractional derivative in Riemann Liouville  fractional derivative sense.

### Example

```julia-repl
julia> fracdiff(x->x^5, 0.5, 2.5, 0.0001, RLDiff_Approx())
141.59707906952633
```

!!! warning
    The RLDiff_Approx algorithm only support for 0 < α < 1.
"""
struct RLDiff_Approx <: RLDiff end


"""
# Riemann Liouville sense derivative using Triangular Strip Matrix to discrete and compute.

    fracdiff(f, α, end_point, h, RLDiff_Matrix())

Using [Triangular Strip Matrix](https://en.wikipedia.org/wiki/Triangle_strip) to approximate fractional derivative.

### Example

```julia-repl
julia> fracdiff(x->x^5, 0.5, 2.5, 0.0001, RLInt_Matrix())
```

!!! info
    Triangular Strip Matrix method returns the derivative in the interval ``[0, T]`` in ```Vector```

```tex
@article{2009,
title={Matrix approach to discrete fractional calculus II: Partial fractional differential equations},
DOI={10.1016/j.jcp.2009.01.014},
author={Podlubny, Igor and Chechkin, Aleksei and Skovranek, Tomas and Chen, YangQuan and Vinagre Jara, Blas M.},
}
```
"""
struct RLDiff_Matrix <: RLDiff end


"""
# Riemann Liouville sense linear spline interpolation.

    fracdiff(f, α, end_point, h, RL_Linear_Spline_Interp())

Using linear spline interpolation method to approximate the Riemann Liouville fractional derivative.
"""
struct RL_Linear_Spline_Interp <: RLDiff end

"""
# Riemann Liouville sense G1 scheme

    fracdiff(f, α, start_point, end_point, h, RL_G1())

Remove the limit symbol in the definition of Grunwald-Letnikov fractional derivative, thereby leading to a discretization scheme in form of truncated series.

!!! tip
        **RL_G1** also can be used to compute fractional integral~
        ``+\\alpha`` for fractional derivative and ``-\\alpha`` for fractional integral.

```tex
@inproceedings{Guo2015FractionalPD,
  title={Fractional Partial Differential Equations and their Numerical Solutions},
  author={Boling Guo and Xueke Pu and Feng-Hui Huang},
  year={2015}
}
```
"""
struct RL_G1 <: RLDiff end

"""
# Riemann Liouville sense D scheme

        fracdiff(f, α, point, h, RL_D())

```tex
@inproceedings{Guo2015FractionalPD,
  title={Fractional Partial Differential Equations and their Numerical Solutions},
  author={Boling Guo and Xueke Pu and Feng-Hui Huang},
  year={2015}
}
```
"""
struct RL_D <: RLDiff end

################################################################
###                    Type definition done                  ###
################################################################


function fracdiff(f::Union{Number, Function}, α, end_point, h, ::RLDiff_Approx)
    #checks(f, α, 0, end_point)
    typeof(end_point) <: Number ? (end_point == 0 ? 0 : f/sqrt(pi*end_point)) : nothing
    summation = 0
    n = Int64(floor(end_point/h))

    @fastmath @inbounds @simd for i ∈ 0:n-1
        summation += (f(end_point-i*h) - f(end_point-(i+1)*h))*((i+1)^(1-α) - i^(1-α))
    end

    result = ((1-α)*f(0)/n^α+summation)*end_point^(-α)*n^α/gamma(2-α)
    return result
end

function fracdiff(f::Union{Number, Function}, α::Float64, end_point::AbstractArray, h, ::RLDiff_Approx)::Vector
    result = map(x->fracdiff(f, α, x, h, RLDiff_Approx()), end_point)
    return result
end




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



#=
Numerical methods for fractional calculus by Li, Changpin
Page 57

Linear Spline Interpolation
=#
function fracdiff(f::Union{Number, Function}, α, x, h, ::RL_Linear_Spline_Interp)
    typeof(x) <: Number ? (x == 0 ? 0 : f/sqrt(pi*x)) : nothing
    N = Int64(floor(x/h))

    result = 0

    @fastmath @inbounds @simd for k = 0:(N+1)
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

function fracdiff(f::Union{Number, Function}, α::Float64, end_point::AbstractArray, h, ::RL_Linear_Spline_Interp)::Vector
    result = map(x->fracdiff(f, α, x, h, RL_Linear_Spline_Interp()), end_point)
    return result
end

function fracdiff(f::Union{Number, Function}, α, start_point, end_point, h, ::RL_G1)
    typeof(end_point) <: Number ? (end_point == 0 ? 0 : f/sqrt(pi*end_point)) : nothing
    N = Int64(floor((end_point-start_point)/h))

    result = zero(Float64)
    @fastmath @inbounds @simd for j = 0:N-1
        result += gamma(j-α)/gamma(j+1)*f(end_point-j*h)
    end

    return h^(-α)/gamma(-α)*result
end

function fracdiff(f::Union{Number, Function}, α, point, h, ::RL_D)
    typeof(end_point) <: Number ? (end_point == 0 ? 0 : f/sqrt(pi*end_point)) : nothing
    N = Int64(floor(point/h))

    result = zero(Float64)
    @fastmath @inbounds @simd for k = 0:N
        result += ωₖₙ(N, k, α)*f(point-k*h)
    end

    return point^(-α)/gamma(-α)*result

end

function ωₖₙ(n, k, α)
    temp = 0

    if k == 0
        temp = -1
    elseif 1 ≤ k ≤ n-1
        temp = 2*k^(1-α)-(k-1)^(1-α)-(k+1)^(1-α)
    elseif k == n
        temp = (α-1)*n^(-α)-(n-1)^(1-α)+n^(1-α)
    end
    return n^α/(α*(1-α))*temp
end