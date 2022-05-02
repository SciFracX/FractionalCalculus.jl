# Symbolic Fractional Differentiation

With the help of [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl) and [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), FractionalCalculus.jl can compute symbolic fractional differentiation and integral now!

First we need to import FractionalCalculus.jl and SymbolicUtils.jl to our workspace:

```julia-repl
julia> using FractionalCalculus, SymbolicUtils
julia> @syms x
```

Then, use ```semidiff``` function to compute symbolic fractional differentiation:

```julia-repl
julia> semidiff(log(x))
log(4x) / sqrt(πx)

julia> semidiff(exp(x))
exp(x)*SpecialFunctions.erf(sqrt(x)) + 1 / sqrt(πx)

julia> semidiff(erf(x))
exp(-0.5x)*SpecialFunctions.besseli(0, 0.5x)
```