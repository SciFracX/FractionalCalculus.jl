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
"""
#abstract type Lubich <:



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
struct Caputo_Direct_First_Diff_known <: Caputo end

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
@book{oldham_spanier_1984,
title={The fractional calculus: Theory and applications of differentiation and integration to arbitrary order},
author={Oldham, Keith B. and Spanier, Jerome},
year={1984}} 
"""

"""

"""
struct GL_Nomenclature <: GL end

"""
"""
struct GL_Lagrange3Interp <: GL end

"""
Riemann Liouville 
"""
struct RLDiff_Approx <: RLDiff end



"""
"""

"""
"""
struct GL_Finite_Difference <: GL end


################################################################
###                    Type defination done                  ###
################################################################

"""
Check if the format of nargins is correct.
"""
function checks(α, start_point, end_point)
    α % 1 != 0 ? nothing : error("α must be a decimal number")
    if typeof(end_point) <: Number
        start_point < end_point ? nothing : DomainError("Domain error! Start point must smaller than end point")
    elseif typeof(end_point) <: AbstractArray
        start_point < minimum(end_point) ? nothing : DomainError("Vector domain error! Start point must smaller than end point")
    else
        ErrorException("Please input correct point you want to compute")
    end
end


"""
# FracDiffAlg

    fracdiff(f::Function, α, start_point, end_point, FracDiffAlg())

### Example

```julia-repl
julia> fracdiff(Function, Order, Start_Point, End_Point, AlgType)
```

Compute the α-order fractional derivative of f from start point to end point with specific algorithm.
"""
function fracdiff(f, α, start_point, end_point, ::FracDiffAlg)
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
# Grünwald–Letnikov sense fractional dervivative.

    fracdiff(f, α, start_point, end_point, GL_Direct())

### Example:

```julia-repl
julia> fracdiff(x->x^5, 0, 0.5, GL_Direct())
```

> Please note Grunwald-Letnikov sense fracdiff only support 0 < α < 1.

Please refer to [Grünwald–Letnikov derivative](https://en.wikipedia.org/wiki/Gr%C3%BCnwald%E2%80%93Letnikov_derivative) for more details.
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
# Caputo sense fractional derivative with first derivative known.

    fracdiff(fd, α, start_point, end_point, Caputo_Direct_First_Diff_known())

If the first order derivative of a function is already known, we can use this method to compute the fractional order derivative more precisely.

The inout function should be the first order derivative of a function

### Example

```julia-repl
julia> fracdiff(x->5*x^4, 0.5, 0, 2.5, Caputo_Direct_First_Diff_known())
```
Return the semi-derivative of ``f(x)=x^5`` at ``x=2.5``.

Compared with **Caputo_Direct** method, this method don't need to specify step size, more precision are guaranteed.
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
# Caputo sense Piecewise algorithm

    fracdiff(f, α, end_point, step_size, Caputo_Piecewise())

Using the **piecewise algorithm** to obtain the fractional derivative at a specific point.

### Example

```julia-repl
julia> fracdiff(x->x^5, 0.5, 2.5, Caputo_Piecewise())
```

Return the fractional derivative of ``f(x)=x^5`` at point ``x=2.5``.

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

    summation = 0
    n = end_point/step_size

    for i in range(0, n, step=1)
        summation +=W₅(i, n, m, α)*first_order(f, i*step_size, step_size)
    end

    result=summation*step_size^(m-α)/gamma(m-α+2)
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

#TODO: Use the improved alg!! This algorithm is not accurate
#This algorithm is not good, still more to do
function fracdiff(f, α, end_point, step_size, ::GL_Nomenclature)
    
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
            append!(ResultArray, fracdiff(f, α, value, step_size, GL_Nomenclature()))
        end
        return ResultArray
    end

    summation = 0
    n = end_point/step_size

    for i in range(0, n-1, step=1)
        summation+=gamma(i-α)/gamma(i+1)*f(end_point-i*end_point/n)
    end
    result=summation*end_point^(-α)*n^α/gamma(-α)
    return result
end

#TODO: This algorithm is same with the above one, not accurate!!!
#This algorithm is not good, still more to do
function fracdiff(f, α, end_point, step_size, ::GL_Lagrange3Interp)
        
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
            append!(ResultArray, fracdiff(f, α, value, step_size, GL_Lagrange3Interp()))
        end
        return ResultArray
    end

    n = end_point/step_size
    summation=0

    for i in range(0, n-1, step=1)
        summation += gamma(i-α)/gamma(i+1)*(f(end_point-i*end_point/n)+1/4*α*(f(end_point-(i-1)*end_point/n)-f(end_point-(i+1)*end_point/n))+1/8*α^2*(f(end_point-(i-1)*end_point/n)-2*f(end_point-i*end_point/n)+f(end_point-(i+1)*end_point/n)))
    end

    result = summation*end_point^(-α)*n^α/gamma(-α)
    return result
end

"""
# Riemann Liouville sense derivative approximation

    fracdiff(f, α, end_point, step_size, RLDiff_Approx())

Using approximation to obtain fractional derivative value.

### Example

```julia-repl
julia> fracdiff(x->x^5, 0.5, 2.5, 0.0001, RLDiff_Approx())
```

!!! warning
    The RLDiff_Approx algorithm only support for 0 < α < 1.

"""
function fracdiff(f, α, end_point, step_size, ::RLDiff_Approx)
        
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
            append!(ResultArray, fracdiff(f, α, value, step_size, RLDiff_Approx()))
        end
        return ResultArray
    end

    summation = 0
    n = end_point/step_size

    for i in range(0, n-1, step=1)
        summation+=(f(end_point-i*end_point/n)-f(end_point-(i+1)*end_point/n))*((i+1)^(1-α)-i^(1-α))
    end

    result=((1-α)*f(0)/n^α+summation)*end_point^(-α)*n^α/gamma(2-α)
    return result
end

function fracdiff(f, α, end_point, step_size, ::GL_Finite_Difference)


    n = end_point/step_size
    result=0

    for i in range(0, n, step=1)
        result+=(-1)^i*gamma(α+1)/(gamma(i+1)*gamma(α-i+1))*f(end_point-i*step_size)
    end

    result1=result/step_size^α
    return result1
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
macro semifracdiff(f, point)
    return :(fracdiff($f, 0.5, $point, 0.0001, RLDiff_Approx()))
end