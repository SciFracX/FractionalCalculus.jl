"""
Riemann-Liouville sense fractional derivative algorithms, please refer to [Riemann-Liouville derivative](https://en.wikipedia.org/wiki/Fractional_calculus#Riemann%E2%80%93Liouville_fractional_derivative)
"""
abstract type RLDiff <: FracDiffAlg end

#RLDiff_LinearInterp maybe??
"""
# Riemann Liouville sense derivative approximation

    fracdiff(f, α, end_point, h, RLDiffL1())

Using the L1 algorithm Linear interpolation to approximate fractional derivative in Riemann Liouville  fractional derivative sense.

### Example

```julia-repl
julia> fracdiff(x->x^5, 0.5, 2.5, 0.0001, RLDiffL1())
141.59707906952633
```

!!! warning
    The RLDiffL1 algorithm only support for 0 < α < 1.
"""
struct RLDiffL1 <: RLDiff end


"""
# Riemann Liouville sense derivative using Triangular Strip Matrix to discrete and compute.

    fracdiff(f, α, end_point, h, RLDiffMatrix())

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
struct RLDiffMatrix <: RLDiff end


"""
# Riemann Liouville sense linear spline interpolation.

    fracdiff(f, α, end_point, h, RLLinearSplineInterp())

Using linear spline interpolation method to approximate the Riemann Liouville fractional derivative.
"""
struct RLLinearSplineInterp <: RLDiff end

"""
# Riemann Liouville sense G1 scheme

    fracdiff(f, α, start_point, end_point, h, RLG1())

Remove the limit symbol in the definition of Grunwald-Letnikov fractional derivative, thereby leading to a discretization scheme in form of truncated series.

!!! tip
        **RLG1** also can be used to compute fractional integral~
        ``+\\alpha`` for fractional derivative and ``-\\alpha`` for fractional integral.

```tex
@inproceedings{Guo2015FractionalPD,
  title={Fractional Partial Differential Equations and their Numerical Solutions},
  author={Boling Guo and Xueke Pu and Feng-Hui Huang},
  year={2015}
}
```
"""
struct RLG1 <: RLDiff end

"""
# Riemann Liouville sense D scheme

        fracdiff(f, α, point, h, RLD())

```tex
@inproceedings{Guo2015FractionalPD,
  title={Fractional Partial Differential Equations and their Numerical Solutions},
  author={Boling Guo and Xueke Pu and Feng-Hui Huang},
  year={2015}
}
```
"""
struct RLD <: RLDiff end

"""
RLDiffL2
"""
struct RLDiffL2 <: RLDiff end

"""
RLDiffL2C
"""
struct RLDiffL2C <: RLDiff end
################################################################
###                    Type definition done                  ###
################################################################

# Riemann Liouville sense L1 method
function fracdiff(f::FunctionAndNumber,
                  α::T,
                  end_point::Real,
                  h::T,
                  ::RLDiffL1) where {T <: Real}
    #checks(f, α, 0, end_point)
    typeof(f) <: Real ? (end_point == 0 ? (return 0) : (return f/sqrt(pi*end_point))) : nothing
    end_point == 0 ? (return 0) : nothing

    summation = zero(T)
    n = floor(Int, end_point/h)

    @fastmath @inbounds @simd for i ∈ 0:n-1
        summation += (f(end_point-i*h) - f(end_point-(i+1)*h))*((i+1)^(1-α) - i^(1-α))
    end

    result = ((1-α)*f(0)/n^α+summation)*end_point^(-α)*n^α/gamma(2-α)
    return result
end

function fracdiff(f::FunctionAndNumber, α::T, end_point::AbstractArray, h::T, ::RLDiffL1)::Vector where {T <: Real}
    result = map(x->fracdiff(f, α, x, h, RLDiffL1()), end_point)
    return result
end




function fracdiff(f::FunctionAndNumber,
                  α::T,
                  end_point::Real,
                  h::T,
                  ::RLDiffMatrix) where {T <: Real}
    N = round(Int, end_point/h+1)
    @views tspan = collect(0:h:end_point)
    return B(N, α, h)*f.(tspan)
end

#Compute the eliminator matrix Sₖ by omiting n-th row
function eliminator(n, row)
    temp = zeros(n, n) + I
    return @views temp[Not(row), :]
end

function B(N, p)
    result = zeros(N, N)
    temp = omega(N, p)

    @inbounds @simd for i ∈ 1:N
        @views result[i, 1:i]=reverse(temp[1:i])
    end

    return result
end
# Multiple dispatch for assigning step size *h*.
function B(N, p, h::T) where {T <: Real}
    result = B(N, p)

    return h^(-p)*result
end




#=
Numerical methods for fractional calculus by Li, Changpin
Page 57

Linear Spline Interpolation
=#
function fracdiff(f::FunctionAndNumber,
                  α::T,
                  x::Real,
                  h::T,
                  ::RLLinearSplineInterp) where {T <: Real}
    typeof(f) <: Real ? (x == 0 ? (return 0) : (return f/sqrt(pi*x))) : nothing
    x == 0 ? (return 0) : nothing
    N = round(Int, x/h)
    result = zero(T)
    @fastmath @inbounds @simd for k = 0:(N+1)
        result += z̄ₘₖ(N, k, α)*f(k*h)
    end
    return 1/(gamma(4-α)h^α)*result
end

function z̄ₘₖ(m::P, k::P, α::T) where {P <: Integer, T <: Real}
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

function fracdiff(f::FunctionAndNumber, α::T, end_point::AbstractArray, h::T, ::RLLinearSplineInterp)::Vector where {T <: Real}
    result = map(x->fracdiff(f, α, x, h, RLLinearSplineInterp()), end_point)
    return result
end

function fracdiff(f::FunctionAndNumber,
                  α::T,
                  start_point::Real,
                  end_point::Real,
                  h::T,
                  ::RLG1) where {T <: Real}
    typeof(f) <: Number ? (end_point == 0 ? (return 0) : (return f/sqrt(pi*end_point))) : nothing
    end_point == 0 ? (return 0) : nothing

    N = round(Int, (end_point-start_point)/h)

    result = zero(T)
    @fastmath @inbounds @simd for j = 0:N-1
        result += gamma(j-α)/gamma(j+1)*f(end_point-j*h)
    end

    return h^(-α)/gamma(-α)*result
end

fracdiff(f::FunctionAndNumber, α, end_point, h::T, ::RLG1) where {T <: Real} = fracdiff(f::FunctionAndNumber, α, 0, end_point, h::Float64, RLG1())

function fracdiff(f::FunctionAndNumber,
                  α::T,
                  point::Real,
                  h::T,
                  ::RLD) where {T <: Real}
    typeof(f) <: Real ? (point == 0 ? (return 0) : (return f/sqrt(pi*point))) : nothing
    point == 0 ? (return 0) : nothing

    N = round(Int, point/h)

    result = zero(T)
    @fastmath @inbounds @simd for k = 0:N
        result += ωₖₙ(N, k, α)*f(point-k*h)
    end

    return point^(-α)/gamma(-α)*result

end

function ωₖₙ(n::P, k::P, α) where {P <: Integer}
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

function fracdiff(f::FunctionAndNumber,
                  α::T,
                  point::Real,
                  h::T,
                  ::RLDiffL2C) where {T <: Real}
    N = round(Int, point/h)
    temp = zero(T)
    for k=-1:N+1
        temp += Ŵ(k, α, N)*f((N-k)*h)
    end
end

function Ŵ(k::P, α::T, N::P) where {P <: Integer, T <: Real}
    expo = 2-α
    if k == -1
        return 1
    elseif k == 0
        return 2^expo-2
    elseif k == 1
        return 3^expo-2^expo
    elseif 2 ≤ k ≤ N-2
        return (k+2)^expo - 2*(k+1)^expo + 2*(k-1)^expo - (k-2)^expo
    elseif k == N-1
        return -N^expo-(N-3)^expo + 2*(N-2)^expo
    elseif k == N
        return -N^expo + 2*(N-1)^expo-(N-2)^expo
    elseif k == N+1
        return N^expo - (N-1)^expo
    end
end
#=
function fracdiff(f::FunctionAndNumber, α, point, h, ::RLDiffL2C)
    
end
=#


function fracdiff(f::FunctionAndNumber,
                  alpha::T,
                  point::Real,
                  h::T,
                  ::RLDiffL2) where {T <: Real}

    n = round(Int, point/h)
    result = zero(T)
    for k=-1:n
        result += WK(k, alpha, n)*f((n-k)*h)
    end

    return f(0)*point^(-alpha)/gamma(1-alpha) + derivative(f, 0)*point^(1-alpha)/gamma(2-alpha) + h^(-alpha)/gamma(3-alpha)*result
end

function WK(k::P, alpha::T, n::P) where {P <: Integer, T <: Real}
    if k == -1
        return 1
    elseif k == 0
        return 2^(2-alpha)-3
    elseif 1 ≤ k ≤ n-2
        return (k+2)^(2-alpha) - 3*(k+1)^(2-alpha) + 3*k^(2-alpha) - (k-1)^(2-alpha)
    elseif k == n-1
        return -2*n^(2-alpha) + 3*(n-1)^(2-alpha) - (n-2)^(2-alpha)
    elseif k == n
        return n^(2-alpha) - (n-1)^(2-alpha)
    end
end