# Fractional integral

```@contents
Pages = ["integral.md"]
```

Fractional integral is defined as follow:

```math
_aD_t^{-\alpha}f(t)=\frac{1}{\Gamma(\alpha)}\int_a^t(t-\tau)^{\alpha-1}f(\tau)d\tau
```
In FractionalCalculus, you can compute the integral of a function with order $\alpha$:

```julia-repl
julia> fracint(x->x, 0.5, 0, 1, 0.0001, RL_Direct())
```

A tuple contains result and estimating error will return

```julia
julia> fracint(x->x, 0.5, 0, 1, 0.0001, RL_())
(0.7522527785271369, 8.022170098417246e-9)
```

## Linear interpolation to approximate function

## Highlight on some algorithms

FractionalCalculus.jl support many algorithms to calculate fractional integral, here, I want to highlight the **Triangular Strip Matrix** method proposed by [Prof Igor](http://people.tuke.sk/igor.podlubny/index.html) to discrete fractional derivative.

```julia-repl
julia> fracdiff(f, α, end_point, h, RLInt_Matrix())
```

With this advancing algorithms, we can not only compute the fractional integral, but also the integer integral! Or more precisely, arbitrary order!!!!🙌 Higher order integral is also a piece of cake!!!!!!🎉

Try to set $\alpha$ as an integer, arbitrary integer of course! I promise you would enjoy it😏