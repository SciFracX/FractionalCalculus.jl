# Fractional integral

```@contents
Pages = ["integral.md"]
```

## Riemann Liouville sense integral

```math
_aD_t^{-\alpha}f(t)=\frac{1}{\Gamma(\alpha)}\int_a^t(t-\tau)^{\alpha-1}f(\tau)d\tau
```
In FractionalCalculus, you can compute the integral of a function with order $\alpha$:

```julia-repl
julia> fracint(x->x, 0.5, 0, 1, 0.0001, RLDirect())
```

A tuple contains result and estimating error will be returned.

```julia
julia> fracint(x->x, 0.5, 0, 1, 0.0001, RLDirect())
(0.7522527785271369, 8.022170098417246e-9)
```