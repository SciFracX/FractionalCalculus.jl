# Fractional integral

```@contents
Pages = ["integral.md"]
```

Fractional integral is defined as follow:

```math
_aD_t^{-\alpha}f(t)=\frac{1}{\Gamma(\alpha)}\int_a^t(t-\tau)^{\alpha-1}f(\tau)d\tau
```
In FractionalCalculus, you can compute the integral of a function with order $\alpha$:

```julia
fracint(x->x, 0.5, 0, 1, 0.0001, RL_Direct())
```

A tuple contains result and estimating error will return

```julia
julia> fracint(x->x, 0.5, 0, 1, 0.0001, RL_())
(0.7522527785271369, 8.022170098417246e-9)
```

## Linear interpolation to approximate function





