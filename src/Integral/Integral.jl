"""
The base type of a fractional integral sense, in FractionalCalculus.jl, all of the fractional integral senses belong to ```FracIntAlg```
"""
abstract type FracIntAlg end


include("RL.jl")

"""
    fracint(f, α, point, h, FracIntAlg())

The general API used for fractional integral computing.
"""
function fracint(f, α, point, h, ::FracIntAlg)
end



## Macros for convenient computing.
"""
    @fracint(f, α, point)

Return the α-order integral of **f** at specific point.

```julia-repl
julia> @fracint(x->x, 0.5, 1)
0.7522525439593486
```
"""
macro fracint(f::Union{Function, Number}, α::Float64, point)
    return :(fracint($f, $α, $point, 0.0001, RLInt_Approx()))
end

"""
    @semifracint(f, point)

Return the semi-integral of **f** at specific point.

```julia-repl
julia> @semifracint(x->x, 1)
0.7522525439593486
```
"""
macro semifracint(f::Union{Function, Number}, point)
    return :(fracint($f, 0.5, $point, 0.0001, RLInt_Approx()))
end