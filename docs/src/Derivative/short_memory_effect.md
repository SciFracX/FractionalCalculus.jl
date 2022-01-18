# Short Memory Effect

In fractional derivative(also in the fractional integral), we already know the derivative of a fixed point depends on its history which is also called the **Non-Local property**. In Grunwald Letnikov sense, when we want to calculate the fractional derivative at a fixed point far away from its starting point, we can neglect its faraway **lower terminal**, and focus on its **recent history** instead.

Which means in the interval $[t-L, t]$, $L$ is the "memory length".

```math
_aD^\alpha_t f(t)\approx _{t-L}D^\alpha_t f(t)
```

By deploying the **Short Memory Effect**, we can reduce our numerical cost while retain the precision in a way.

!!! info "Short memory effect in FDE"
    Want to see how the short memory is used in fractional differential equations? Please see [short memory effect in FDE](http://scifracx.org/FractionalDiffEq.jl/dev/system_of_FDE/#Short-memory-effect-in-FDE)