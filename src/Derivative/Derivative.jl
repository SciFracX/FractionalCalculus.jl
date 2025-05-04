"""
Base type of fractional differentiation algorithms, in FractionalCalculus.jl, all of the fractional derivative algorithms belong to ```FracDiffAlg```
"""
abstract type FracDiffAlg end

const FunctionAndNumber = Union{Function, Number}

include("Caputo.jl")
include("GL.jl")
include("Hadamard.jl")
include("Riesz.jl")
include("RL.jl")
include("CaputoFabrizio.jl")
include("ABC.jl")
include("Chebyshev.jl")


"""
    fracdiff(f, α, point, h, FracDiffAlg())

```fracdiff``` is the general function for computing fractional derivative.

### Example

```julia-repl
julia> fracdiff(x->x^2, 0.5, 1, 0.001, RLDiffL1())
1.5044908143658473
```

There are many algorithms can be used to computing fractional derivative, please refer to our manual for more details: http://scifracx.org/FractionalCalculus.jl/dev/Derivative/derivativeapi/.
"""
function fracdiff(f, α, point, h, ::FracDiffAlg)

end

function checkfracdiff(f, α, point, h)
    point > 0 ? nothing : DomainError("Only positive points")
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
    return :(fracdiff($f, $α, $point, 0.0001, RLDiffL1()))
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
    return :(fracdiff($f, 0.5, $point, 0.0001, RLDiffL1()))
end

