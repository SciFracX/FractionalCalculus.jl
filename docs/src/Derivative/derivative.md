# Fractional derivative

```@contents
Pages = ["derivative.md"]
```

To get start with fractional derivative, you need to know that unlike Newtonian derivatives, fractional derivative is defined via integral.

!!! tip "Non-Local Operators"
	It is noteworthy that the fractional derivatives are not local operators, which means that we cannot calculate the fractional derivative solely on the basis of function values of $f(x)$ taken from neighborhood of the point $x$ where er work. Instead, we have to take the entire history of $f(x)$ (i.e., all function values f(x) for $0<x<a$) into account.

## Riemann Liouville sense derivative

Riemann Liouville sense derivative is bult upon the Riemann Liouville sense integral.
```math
_aD^\alpha_tf(t)=\frac{d^n}{dt^n}\ _aD^{-(n-\alpha)}_tf(t)=\frac{d^n}{dt^n}\ _aI^{n-\alpha}_tf(t)
```

```math
_tD^\alpha_bf(t)=\frac{d^n}{dt^n}\ _tD^{-(n-\alpha)}_bf(t)=\frac{d^n}{dt^n}\ _tI^{n-\alpha}_bf(t)
```

We can use FractionalCalculus.jl to compute fractional derivative:

```julia
fracdiff(x->x, 0.5, 1, 0.0001, RLDiff_Approx())
```

## Caputo sense derivative

In **FractionalCalculus.jl**, let's see, if you want to calculate the $0.4$ order fractional derivative of $f(x)=x$ at a $x=1$ with step size $0.0001$, simply typing these:


```julia
fracdiff(x->x, 0.4, 1, 0.0001, Caputo_Piecewise())
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

## Grünwald Letnikov sense derivative

```math
D^\alpha f(t)=\displaystyle \lim_{h\rightarrow0}\frac{1}{h^\alpha}\sum_{0\leq m\lt\infty}(-1)^m {{\alpha}\choose{m}}f(t+(\alpha-m)h)
```