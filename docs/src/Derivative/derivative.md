# Fractional derivative

To get start with fractional derivative, you need to know the fractional derivative is defined via integral
```math
D^{\alpha}_{t}f(t)=\frac{1}{\Gamma(n-\alpha)}\int^t_0
```


```julia
fracdiff(x->x, 0.5, 0, 1, Caputo())
```

There many types of definitions of fractional derivative, Caputo is one of these useful definitions. The Caputo fractional derivative is first be proposed in [Michele Caputo's Paper](https://doi.org/10.1111/j.1365-246X.1967.tb02303.x), 
```math
^CD_t^\alpha f(t) = \frac{1}{\Gamma(n-\alpha)}\int_0^t\frac{f^{(n)}(\tau)d\tau}{(t-\tau)^{\alpha+1-n}}, n=\lceil{\alpha}\rceil
```

In FractionalCalculus.jl, we use **complex step differentiation**, which is also called **CSD**, to obtain the numerical solution in $0<\alpha<1$, so here the value of n is $1$.

> We use **Cleve Moler's** annotation about **complex step differentiation**:
>
> > **Complex step differentiation** is a technique that employs complex arithmetic to obtain the numerical value of the first derivative of a real valued analytic function of a real variable, avoiding the loss of precision inherent in traditional finite differences.

We use the above information to rewrite the equation:
```math
^CD^{\alpha}_tf(t) = \frac{1}{\Gamma(1-\alpha)}\int_0^t\frac{f'(\tau)d\tau}{(t-\tau)^\alpha}=\frac{1}{\Gamma(1-\alpha)}\int_0^t(t-\tau)^{-\alpha}\frac{f(\tau+h)-f(\tau)}{h} d\tau
```
By applying **complex step differentiation**, $f'(\tau)\approx\frac{Im\{f(\tau+ih)\}}{h}$, now we get:
```math
^CD^\alpha_tf(t)=\frac{1}{\Gamma(1-\alpha)}\int^t_0(t-\tau)^{-\alpha}\frac{Im\{f(\tau+ih)\}}{h}d\tau
```

