"""
Caputo-Fabrizio sense fractional derivative algorithms
"""
abstract type CaputoFabrizio <: FracDiffAlg end

"""
# Atangana Seda algorithm for approximating Caputo-Fabrizio fractional Derivatives

    fracdiff(f, α, point, h, CaputoFabrizioAS())

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.001, CaputoFabrizioAS())
0.9893315392351845
```

### Reference

```tex
Fractional Stochastic Differential Equations
```
"""
struct CaputoFabrizioAS <: CaputoFabrizio end

function fracdiff(f, α, point, h, ::CaputoFabrizioAS)
    M = 1-α+α/gamma(α)
    n::Int = round(Int, point/h)
    result = zero(Float64)
    for j in 0:n
        result += (f((j+1)*h)-f(j*h))/h*((1-α)/α*(exp(-α/(1-α)*(n-j)*h)-exp(-α/(1-α)*(n-j+1)*h)))
    end
    return result*M/(1-α)
end