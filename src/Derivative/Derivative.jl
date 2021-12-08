# The structure of the fractional derivative:
# (1) Type definition.
# (2) Checks function to ensure the input parameters are correct and handle some special cases, for example, 0 end point case.
# (3) Multiple dispatch for fracdiff of different algorithms, the first is used to compute the value of a specific point, the second is used to compute the fractional derivative of a Vector, other functions are auxillary functions to help with the coefficients.
# (4) Fractional derivative macros for convenient computing.

using QuadGK, SpecialFunctions, SpecialMatrices
using LinearAlgebra, InvertedIndices

"""
Base type of fractional differentiation algorithms, in FractionalCalculus.jl, all of the fractional derivative algorithms belong to ```FracDiffAlg```
"""
abstract type FracDiffAlg end

"""
Caputo sense fractional derivative algorithms, please refer to [Caputo derivative](https://en.wikipedia.org/wiki/Fractional_calculus#Caputo_fractional_derivative) for more details.
"""
abstract type Caputo <: FracDiffAlg end

"""
Grünwald–Letnikov sense fractional derivative algorithms, please refer to [Grünwald–Letnikov derivative](https://en.wikipedia.org/wiki/Gr%C3%BCnwald%E2%80%93Letnikov_derivative) for more details
"""
abstract type GL <: FracDiffAlg end

"""
Riemann-Liouville sense fractional derivative algorithms, please refer to [Riemann-Liouville derivative](https://en.wikipedia.org/wiki/Fractional_calculus#Riemann%E2%80%93Liouville_fractional_derivative)
"""
abstract type RLDiff <: FracDiffAlg end

"""
Hadamard sense fractional derivative algorithms.
"""
abstract type Hadamard <: FracDiffAlg end

"""
Riesz sense fractional derivative algorithms.
"""
abstract type Riesz <: FracDiffAlg end




"""
Note *Caputo Direct* algorithms belong to direct computing, precise are ensured, but maybe cause more memory allocation and take more compilation time.

Using the direct mathematic expression:

```math
^CD_t^α=\\frac{1}{\\Gamma(n-α)}\\int_0^t\\frac{f^{(n)}(τ)}{(t-τ)^{α+1-n}}
```

As for the derivative inside the integral, we use the **Complex Step Differentiation** to obtain a more accurate value.
"""
struct Caputo_Direct <: Caputo end

"""
With first derivative known, we can direct use the derivative to compute a more accurate result.
"""
struct Caputo_Direct_First_Diff_Known <: Caputo end

"""
With first and second derivative known, we can direct use derivative and second derivative to compute a more accurate result.
"""
struct Caputo_Direct_First_Second_Diff_Known <: Caputo end

"""
Grunwald Letnikov direct compute method to obtain fractional derivative, precision are guaranteed but cause more memory allocation and compilation time.
"""
struct GL_Direct <: GL end


"""
@article{LI20113352,
title = {Numerical approaches to fractional calculus and fractional ordinary differential equation},
author = {Changpin Li and An Chen and Junjie Ye},
}
"""

"""
Using piecewise linear interpolation function to approximate input function and combining Caputo derivative then implement summation.
"""
struct Caputo_Piecewise <: Caputo end


"""
@article{10.1093/imanum/17.3.479,
author = {DIETTELM, KAI},
title = "{Generalized compound quadrature formulae for finite-part integrals}",
journal = {IMA Journal of Numerical Analysis},
doi = {10.1093/imanum/17.3.479},
}
"""

"""
Diethelm's method for computingn Caputo sense fractional derivative using [product trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule) and [Richardsin extrapolation](https://en.wikipedia.org/wiki/Richardson_extrapolation)
"""
struct Caputo_Diethelm <: Caputo end


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
Using **Grünwald–Letnikov finite difference method** to compute Grünwald–Letnikov sense fractional derivative.
"""
struct GL_Finite_Difference <: GL end


"""
@article{article,
author = {Yin, Chuntao and Li, Changpin and Bi, Qinsheng},
year = {2018},
month = {01},
pages = {},
title = {Approximating Hadamard Derivative and Fractional Differential Equation via the Finite Part Integral},
journal = {SSRN Electronic Journal},
doi = {10.2139/ssrn.3281675}
}
"""

"""
Using finite part integral of **left rectangular formula** to compute the fractional hadamard derivative.
"""
struct Hadamard_LRect <: Hadamard end

"""
Using finite part integral of **right rectangular formula** to compute the fractional hadamard derivative.
"""
struct Hadamard_RRect <: Hadamard end

"""
Using finite part integral of **trapezoidal formula** to compute the fractional hadamard derivative.
"""
struct Hadamard_Trap <: Hadamard end


"""
Using the Triangular Strip Matrices to compute the fractional derivative Riesz derivative.
"""
struct Riesz_Symmetric <: Riesz end

################################################################
###                    Type defination done                  ###
################################################################

"""
Check if the format of nargins is correct.
"""
function checks(f::Union{Function, Number}, α::Float64, start_point::Float64, end_point::Float64)
    α % 1 != 0 ? nothing : error("α must be a decimal number")
    if isa(end_point, Number)
        start_point < end_point ? nothing : DomainError("Domain error! Start point must smaller than end point")
    elseif isa(end_point, AbstractArray)
        start_point < minimum(end_point) ? nothing : DomainError("Vector domain error! Start point must smaller than end point")
    else
        ErrorException("Please input correct point you want to compute")
    end

    end_point == 0 ? 0 : nothing
end

function checks(f::Number, α::Float64, start_point::Float64, end_point::Float64)
    f = f == 0 ? 0 : f/sqrt(pi*end_point)
end





"""
# Caputo sense fractional derivative.

    fracdiff(f::Function, α, start_point, end_point, step_size ::Caputo)

### Example: 

```julia-repl
julia> fracdiff(x->x^5, 0.5, 0, 2.5, 0.0001, Caputo_Direct())
```

Returns a tuple (**result**, **error**), which means the value of this derivative is 141.59714979764541, and the error estimate is 1.1532243848672914e-6.

When the input end points is an array, **fracdiff** will compute 

Refer to [Caputo derivative](https://en.wikipedia.org/wiki/Fractional_calculus#Caputo_fractional_derivative)
"""
function fracdiff(f::Union{Function, Number}, α::Float64, start_point, end_point, h::Float64, ::Caputo_Direct)
    checks(f, α, start_point, end_point)

    #Using complex step differentiation to calculate the first order differentiation
    g(τ) = imag(f(τ+1*im*h) ./ h) ./ ((end_point-τ) .^α)
    result = quadgk(g, start_point, end_point) ./gamma(1-α)
    return result
end

function fracdiff(f::Union{Number, Function}, α::Float64, start_point, end_point::AbstractArray, h::Float64, ::Caputo_Direct)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, start_point, value, h, Caputo_Direct()))
    end
    return ResultArray
end



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
    checks(f, α, start_point, end_point)

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



"""
# Caputo sense fractional derivative with first derivative known.

    fracdiff(fd, α, start_point, end_point, Caputo_Direct_First_Diff_Known())

If the first order derivative of a function is already known, we can use this method to compute the fractional order derivative more precisely.

The inout function should be the first order derivative of a function

### Example

```julia-repl
julia> fracdiff(x->5*x^4, 0.5, 0, 2.5, Caputo_Direct_First_Diff_Known())
```
Return the semi-derivative of ``f(x)=x^5`` at ``x=2.5``.

Compared with **Caputo_Direct** method, this method don't need to specify step size, more precision are guaranteed.
"""
function fracdiff(fd::Function, α, start_point, end_point, ::Caputo_Direct_First_Diff_Known)

    g(τ) = (end_point-τ) .^ (-α) * fd(τ)
    result = quadgk(g, start_point, end_point) ./gamma(1-α)    
    return result
end

function fracdiff(f::Function, α, start_point, end_point::AbstractArray, step_size, ::Caputo_Direct_First_Diff_Known)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, start_point, value, step_size, Caputo_Direct_First_Diff_Known()))
    end
    return ResultArray
end



"""
# Caputo sense fractional derivative with first and second derivative known.

    fracdiff(fd1, fd2, α, start_point, end_point, Caputo_Direct_First_Second_Diff_Known)

If the first and second order derivative of a function is already known, we can use this method to compute the fractional order derivative more precisely.


### Example

```julia-repl
julia> fracdiff(x->5*x^4, x->20*x^3, 0.5, 0, 2.5, Caputo_Direct_First_Second_Diff_known())
```
Return the semi-derivative of ``f(x)=x^5`` at ``x=2.5``.
    
Compared with **Caputo_Direct** method, this method don't need to specify step size, more precision are guaranteed.
"""
function fracdiff(fd1::Function, fd2, α, start_point, end_point, ::Caputo_Direct_First_Second_Diff_Known)

    temp1 = fd1(start_point) .* (end_point-start_point) .^ (1-α) ./ gamma(2-α)
    g(τ) = fd2(τ) .* ((end_point-τ) .^ (1-α))
    temp2 = quadgk(g, start_point, end_point) ./ gamma(2-α)
    result = temp1 .+ temp2
    return result
end

function fracdiff(f::Function, α::Float64, start_point, end_point::AbstractArray, h, ::Caputo_Direct_First_Second_Diff_Known)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, start_point, value, h, Caputo_Direct_First_Second_Diff_Known()))
    end
    return ResultArray
end


"""
# Caputo sense Piecewise algorithm

    fracdiff(f, α, end_point, h, Caputo_Piecewise())

Using the **piecewise algorithm** to obtain the fractional derivative at a specific point.

### Example

```julia-repl
julia> fracdiff(x->x^5, 0.5, 2.5, 0.001, Caputo_Piecewise())
```

Return the fractional derivative of ``f(x)=x^5`` at point ``x=2.5``.
"""
function fracdiff(f::Union{Function, Number}, α::Float64, end_point, h, ::Caputo_Piecewise)
    checks(f, α, 0, end_point)

    m = floor(α)+1

    summation = 0
    n = Int64(end_point/h)

    for i in range(0, n, step=1)
        summation +=W₅(i, n, m, α)*first_order(f, i*h, h)
    end

    result=summation*h^(m-α)/gamma(m-α+2)
    return result
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

function fracdiff(f::Union{Number, Function}, α::Float64, end_point::AbstractArray, h, ::Caputo_Piecewise)::Vector
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, value, h, Caputo_Piecewise()))
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
    checks(f, α, 0, end_point)

    summation = 0
    n = Int64(end_point/h)

    for i in range(0, n-1, step=1)
        summation+=gamma(i-α)/gamma(i+1)*f(end_point-i*end_point/n)
    end
    result=summation*end_point^(-α)*n^α/gamma(-α)
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
    checks(f, α, 0, end_point)

    n = end_point/h
    summation=0

    for i in range(0, n-1, step=1)
        summation += gamma(i-α)/gamma(i+1)*(f(end_point-i*end_point/n)+1/4*α*(f(end_point-(i-1)*end_point/n)-f(end_point-(i+1)*end_point/n))+1/8*α^2*(f(end_point-(i-1)*end_point/n)-2*f(end_point-i*end_point/n)+f(end_point-(i+1)*end_point/n)))
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
    n = end_point/h

    @fastmath @inbounds @simd for i in range(0, n-1, step=1)
        summation+=(f(end_point-i*end_point/n)-f(end_point-(i+1)*end_point/n))*((i+1)^(1-α)-i^(1-α))
    end

    result=((1-α)*f(0)/n^α+summation)*end_point^(-α)*n^α/gamma(2-α)
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
    checks(f, α, 0, end_point)

    n = end_point/h
    result=0

    @fastmath @simd for i in range(0, n, step=1)
        result+=(-1)^i*gamma(α+1)/(gamma(i+1)*gamma(α-i+1))*f(end_point-i*h)
    end

    result1=result*h^(-α)
    return result1
end

function fracdiff(f::Union{Function, Number}, α::AbstractArray, end_point, h, ::GL_Finite_Difference)::Vector
    ResultArray = Float64[]

    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracint(f, α, value, h, GL_Finite_Difference()))
    end
    return ResultArray    
end

"""
# Caputo sense Diethelm computation
    
    fracdiff(f, α, end_point, h, Caputo_Diethelm())

Using quadrature weights(derived from product trapezoidal rule) to approximate the derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.007, Caputo_Diethelm())
1.128378318687192
```

!!!info "Scope"
    0 < α < 1
"""
function fracdiff(f::Union{Function, Number}, α::Float64, end_point, h, ::Caputo_Diethelm)
    checks(f, α, 0, end_point)

    N = end_point/h

    result=0

    @fastmath @inbounds @simd for i in range(0, N, step=1)
        result+=quadweights(i, N, α)*(f(end_point-i*h)-f(0))
    end
    return result*h^(-α)/gamma(2-α)

end
function quadweights(n, N, α)
    if n == 0
        return 1
    elseif  0 < n < N
        return (n+1)^(1-α)-2*n^(1-α)+(n-1)^(1-α)
    elseif n == N
        return (1-α)*N^(-α)-N^(1-α)+(N-1)^(1-α)
    end
end
function nderivative(f, n, end_point, h)
    temp = 0

    @fastmath @inbounds @simd for k in range(0, n, step=1)
        temp += (-1)^(k+n)*binomial(n, k)*f(end_point+k*h)
    end
    return h^(-n)*temp
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
    N=Int64(end_point/h+1)
    @views tspan=collect(0:h:end_point)
    return B(N, α, h)*f.(tspan)
end

#Compute the eliminator matrix Sₖ by omiting n-th row
function eliminator(n, row)
    temp = zeros(n, n)+I
    return @views temp[Not(row), :]
end

#Generate elements in Matrix.
function omega(n, p)
    omega = zeros(n+1)

    omega[1]=1
    @fastmath @inbounds @simd for i in range(1, n, step=1)
        omega[i+1]=(1-(p+1)/i)*omega[i]
    end
    
    return omega

end
function B(N, p, h::Float64)
    result=zeros(N, N)
    temp=omega(N, p)

    @inbounds @simd for i in range(1, N, step=1)
        @views result[i, 1:i]=reverse(temp[1:i])
    end

    return h^(-p)*result
end
# Multiple dispatch for not assigning step size *h*.
function B(N, p)
    result=zeros(N, N)
    temp=omega(N, p)

    @inbounds @simd for i in range(1, N, step=1)
        @views result[i, 1:i]=reverse(temp[1:i])
    end

    return result
end


"""
    fracdiff(f, α, a, x, h, Hadamard_LRect())

Compute Hadamard sense using left rectangular formula.
"""
function fracdiff(f, α, x₀, x, h, ::Hadamard_LRect)
    N=Int64((x-x₀)/h)
    result=0

    for i in range(0, N-1, step=1)
        result+=LRectCoeff(i, N, h, α, x₀)*f(x₀+i*h)
    end

    return result
end
function LRectCoeff(i, n, h, α, x₀)
    if 0 ≤ i ≤ n-2
        return 1/gamma(1-α)*((log((x₀+n*h)/(x₀+i*h)))^(-α) - (log((x₀+n*h)/(x₀+(i+1)*h)))^(-α))
    elseif i == n-1
        return 1/gamma(1-α)*(log((x₀+n*h)/(x₀+i*h)))^(-α)
    end
end

"""
    fracdiff(f, α, a, x, h, Hadamard_Trap())

Compute Hadamard sense using right rectangular formula.
"""
function fracdiff(f, α, x₀, x, h, ::Hadamard_RRect)
    N=Int64((x-x₀)/h)
    result=0

    for i in range(0, N-1, step=1)
        result+=LRectCoeff(i, N, h, α, x₀)*f(x₀+(i+1)*h)
    end


    return result
end
function RRectCoeff(i, n, h, α, x₀)
    if i == 0
        return 1/2*LRectCoeff(i, n, h, α, x₀)
    elseif 1 ≤ i ≤ n-1
        return 1/2*(LRectCoeff(i-1, n, h, α, x₀)+LRectCoeff(i, n, h, α, x₀))
    elseif i == n
        return 1/2*LRectCoeff(i-1, n, h, α, x₀)
    end
end

"""
    fracdiff(f, α, a, x, h, Hadamard_Trap())

Compute Hadamard sense using Trapezoidal formula.
"""
function fracdiff(f, α, x₀, x, h, ::Hadamard_Trap)
    N=Int64((x-x₀)/h)

    result=0

    for i in range(0, N, step=1)
        result+=RRectCoeff(i, N, h, α, x₀)*f(x₀+i*h)
    end

    return result
end

"""
    fracdiff(f, α, end_point, h, Riesz_Symmetric())

Compute fractional derivative of Riesz sense using Triangular Strip Matrix algorithm.
"""
function fracdiff(f, α, end_point, h, ::Riesz_Symmetric)
    N=Int(end_point/h)

    mat=RieszMatrix(α, N+1, h)
    return mat*f.(collect(0:h:end_point))
end
function RieszMatrix(α, N, h)
    caputo=B(N+1, α)
    caputo=caputo[2:(N+1), 1:N]
    result=1/2*(caputo+caputo')
    result=h^(-α)*result

    return result
end



"""
P-th precision polynomial generate function

```math
g_p(z)=\\sum_{k=1}^p \\frac{1}{k}(1-z)^k
```
"""
function genfun(p)
    a=collect(1:p+1)
    A=Vandermonde(a)'
    return (1 .-a')*inv(A')
end

"""
"""
function getvec(α, n, g)
    p=length(g)-1
    b=1+α
    g0=g[1]
    w=zeros(1, n)
    w[1]=g[1]^α

    for m=2:p
        M=m-1
        dA=b/M
        w[m]=-collect((g[2:m].*collect(((1-dA):-dA:(1-p*dA)))))*w[M:-1:1]'/g0
    end

    for k=p+1:n
        M=k-1
        dA=b/M
        w[k]=-collect(g[2:(p+1)].*collect((1-dA):-dA:(1-b)))*w[M:-1:1]'/g0
    end
end


## Macros for convenient computing.
"""
    @fracdiff(f, α, point)

Return the α-order derivative of **f** at specific point.

```julia-repl
julia> @fracdiff(x->x, 0.5, 1)
1.1283791670955188
```
"""
macro fracdiff(f, α, point)
    return :(fracdiff($f, $α, $point, 0.0001, RLDiff_Approx()))
end

"""
    @semifracdiff(f, point)

Return the semi-derivative of **f** at spedific point.

```julia-repl
julia> @semifracdiff(x->x, 1)
1.1283791670955188
```
"""
macro semifracdiff(f::Union{Number, Function}, point)
    return :(fracdiff($f, 0.5, $point, 0.0001, RLDiff_Approx()))
end


