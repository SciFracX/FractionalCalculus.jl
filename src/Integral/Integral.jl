"""
The base type of a fractional integral sense, in FractionalCalculus.jl, all of the fractional integral senses belong to ```FracIntAlg```
"""
abstract type FracIntAlg end


include("RL.jl")
include("Hadamard.jl")

include("SymbolicInt.jl")

"""
    fracint(f, α, point, h, FracIntAlg())

The general function used for fractional integral computing.

### Example

```julia-repl
julia> fracint(x->x^2, 0.5, 1, 0.001, RLInt_Approx())
0.6017877646029814
```

There are many algorithms can be used to computing fractional integral, please refer to our manual for more details: [manual](http://scifracx.org/FractionalCalculus.jl/dev/Integral/integralapi/).
"""
function fracint(f, α, point, h, ::FracIntAlg)
end


## Macros for convenient computing.
"""
    @fracint(f, α, point)

Return the α-order integral of ``f`` at specific point.

```julia-repl
julia> @fracint(x->x, 0.5, 1)
0.7522525439593486
```
"""
macro fracint(f, α, point)
    return :(fracint($f, $α, $point, 0.0001, RLIntApprox()))
end

"""
    @semifracint(f, point)

Return the semi-integral of ``f`` at specific point.

```julia-repl
julia> @semifracint(x->x, 1)
0.7522525439593486
```
"""
macro semifracint(f, point)
    return :(fracint($f, 0.5, $point, 0.0001, RLIntApprox()))
end