import FractionalCalculus.FracDiffAlg

"""
Caputo sense fractional derivative algorithms, please refer to [Caputo derivative](https://en.wikipedia.org/wiki/Fractional_calculus#Caputo_fractional_derivative) for more details.
"""
abstract type Caputo <: FracDiffAlg end


"""
# Caputo sense fractional derivative.

    fracdiff(f::Function, α, start_point, end_point, step_size, ::Caputo)

### Example: 

```julia-repl
julia> fracdiff(x->x^5, 0.5, 0, 2.5, 0.0001, Caputo_Direct())
```

Returns a tuple (**result**, **error**), which means the value of this derivative is 141.59714979764541, and the error estimate is 1.1532243848672914e-6.

Refer to [Caputo derivative](https://en.wikipedia.org/wiki/Fractional_calculus#Caputo_fractional_derivative)

!!! note
    ```Caputo Direct``` algorithms belong to direct computing, precise are ensured, but maybe cause more memory allocation and take more compilation time.

Using the direct mathematic expression:

```math
^CD_t^α=\\frac{1}{\\Gamma(n-α)}\\int_0^t\\frac{f^{(n)}(τ)}{(t-τ)^{α+1-n}}
```

As for the derivative inside the integral, we use the **Complex Step Differentiation** to obtain a more accurate value.
"""
struct Caputo_Direct <: Caputo end

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
struct Caputo_Direct_First_Diff_Known <: Caputo end

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
struct Caputo_Direct_First_Second_Diff_Known <: Caputo end


"""
# Caputo sense Piecewise algorithm

    fracdiff(f, α, end_point, h, Caputo_Piecewise())

Using piecewise linear interpolation function to approximate input function and combining Caputo derivative then implement summation.

### Example

```julia-repl
julia> fracdiff(x->x^5, 0.5, 2.5, 0.001, Caputo_Piecewise())
```

Return the fractional derivative of ``f(x)=x^5`` at point ``x=2.5``.

```tex
@article{LI20113352,
title = {Numerical approaches to fractional calculus and fractional ordinary differential equation},
author = {Changpin Li and An Chen and Junjie Ye},
}
```

"""
struct Caputo_Piecewise <: Caputo end



"""
# Caputo sense Diethelm computation
    
    fracdiff(f, α, end_point, h, Caputo_Diethelm())

Using quadrature weights(derived from product trapezoidal rule) to approximate the derivative.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.007, Caputo_Diethelm())
1.128378318687192
```

!!! info "Scope"
    0 < α < 1

Diethelm's method for computingn Caputo sense fractional derivative using [product trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule) and [Richardsin extrapolation](https://en.wikipedia.org/wiki/Richardson_extrapolation)

```tex
@article{10.1093/imanum/17.3.479,
author = {DIETTELM, KAI},
title = "{Generalized compound quadrature formulae for finite-part integrals}",
journal = {IMA Journal of Numerical Analysis},
doi = {10.1093/imanum/17.3.479},
}
```

"""
struct Caputo_Diethelm <: Caputo end


"""
# Caputo sense fractioal derivative with p-th order precision.

    fracdiff(f, α, t, p, Caputo_High_Precision())

Use the high precision algorithm to compute the Caputo sense fractional derivative.

The **p** here is the grade of precision.
"""
struct Caputo_High_Precision <: Caputo end


################################################################
###                    Type definition done                  ###
################################################################



function fracdiff(f::Union{Function, Number}, α::Float64, start_point, end_point, h::Float64, ::Caputo_Direct)
    #checks(f, α, start_point, end_point)
    typeof(f) <: Number ? (end_point == 0 ? 0 : f/sqrt(pi*end_point)) : nothing
    #Using complex step differentiation to calculate the first order differentiation
    g(τ) = imag(f(τ+1*im*h) ./ h) ./ ((end_point-τ) .^α)
    result = quadgk(g, start_point, end_point) ./gamma(1-α)
    return result
end

function fracdiff(f::Union{Number, Function}, α::Float64, start_point, end_point::AbstractArray, h::Float64, ::Caputo_Direct)::Vector
    #=
    ResultArray = Float64[]
    for (_, value) in enumerate(end_point)
        append!(ResultArray, fracdiff(f, α, start_point, value, h, Caputo_Direct()))
    end
    return ResultArray
    =#
    result = map(x->fracdiff(f, α, start_point, x, h, Caputo_Direct())[1], end_point)
    return result
end

#=
These might be redundant, leave it here for now, maybe useful one day ¯\_(ツ)_/¯
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
=#

function fracdiff(f::Union{Function, Number}, α::Float64, end_point, h, ::Caputo_Piecewise)
    typeof(f) <: Number ? (end_point == 0 ? (return 0) : (return f/sqrt(pi*end_point))) : nothing
    end_point == 0 ? (return 0) : nothing
    m = floor(α)+1

    summation = zero(Float64)
    n = Int64(floor(end_point/h))

    @fastmath @inbounds @simd for i ∈ 0:n
        summation += W₅(i, n, m, α)*first_order(f, i*h, h)
    end

    result = summation*h^(m-α)/gamma(m-α+2)
    return result
end
function W₅(i, n, m, α)
    if i == 0
        return n^(m-α)*(m-α+1-n) + (n-1)^(m-α+1)
    elseif i == n
        return 1
    else
        (n-i-1)^(m-α+1) + (n-i+1)^(m-α+1) - 2*(n-i)^(m-α+1)
    end
end
# Deploy Complex Step Differentiation to compute the first order derivative.
function first_order(f, point, h)
    return imag(f(point + im*h))/h
end


function fracdiff(f::Union{Number, Function}, α::Float64, end_point::AbstractArray, h, ::Caputo_Piecewise)::Vector
    result = map(x->fracdiff(f, α, x, h, Caputo_Piecewise()), end_point)
    return result
end


#=
Caputo Diethelm algorithm
=#
function fracdiff(f::Union{Function, Number}, α::Float64, end_point::Number, h, ::Caputo_Diethelm)
    #checks(f, α, 0, end_point)
    typeof(f) <: Number ? (end_point == 0 ? 0 : f/sqrt(pi*end_point)) : nothing
    N = Int64(floor(end_point/h))

    result = zero(Float64)

    @fastmath @inbounds @simd for i ∈ 0:N
        result += quadweights(i, N, α)*(f(end_point-i*h) - f(0))
    end

    return result*h^(-α)/gamma(2-α)
end
function quadweights(n, N, α)
    if n == 0
        return 1
    elseif  0 < n < N
        return (n+1)^(1-α) - 2*n^(1-α) + (n-1)^(1-α)
    elseif n == N
        return (1-α)*N^(-α) - N^(1-α) + (N-1)^(1-α)
    end
end

function fracdiff(f, α::Float64, end_point::AbstractArray, h, ::Caputo_Diethelm)::Vector
    result = map(x->fracdiff(f, α, x, h, Caputo_Diethelm()), end_point)
    return result
end


#=
Caputo sense high precision algorithm
=#
function fracdiff(y, α, t, p, ::Caputo_High_Precision)
    h = t[2]-t[1]
    t = t[:]
    n = length(t)
    y = y.(t)
    q = Int64(ceil(α))
    r = Int64(max(p, q))
    R = reverse(Vandermonde(collect(0:(r-1)).*h))
    c = inv(R)*y[1:r]
    u = 0
    du = 0
    for i = 1:r
        u = u.+c[i]*t.^(i-1)
    end
    if q < r
        for i = (q+1):p
            du = du.+c[i]*t.^(i-1-α)*gamma(1)/gamma(i-α)
        end
    end

    v = y[:] - u[:]
    dv = fracdiff(v, α, t, p, GL_High_Precision())
    dy = dv[:] + du[:]
    return dy
end