# Fractional derivative

```@contents
Pages = ["derivative.md"]
```

To get start with fractional derivative, you need to know that unlike Newtonian derivatives, fractional derivative is defined via integral.

!!! tip "Non-Local Operators"
	It is noteworthy that the fractional derivatives are not local operators, which means that we cannot calculate the fractional derivative solely on the basis of function values of $f(x)$ taken from neighborhood of the point $x$. Instead, we have to take the entire history of $f(x)$ (i.e., all function values $f(x)$ for $0<x<a$) into account.

## Riemann Liouville sense derivative

Riemann Liouville sense derivative is bult upon the Riemann Liouville sense integral.
```math
_aD^\alpha_tf(t)=\frac{d^n}{dt^n}\ _aD^{-(n-\alpha)}_tf(t)=\frac{d^n}{dt^n}\ _aI^{n-\alpha}_tf(t)
```

```math
_tD^\alpha_bf(t)=\frac{d^n}{dt^n}\ _tD^{-(n-\alpha)}_bf(t)=\frac{d^n}{dt^n}\ _tI^{n-\alpha}_bf(t)
```

We can use FractionalCalculus.jl to compute fractional derivative:

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.0001, RLDiff_Approx())
```

## Caputo sense derivative

In **FractionalCalculus.jl**, let's see, if you want to calculate the $0.4$ order fractional derivative of $f(x)=x$ at a $x=1$ with step size $0.0001$, simply typing these:


```julia-repl
julia>fracdiff(x->x, 0.4, 1, 0.0001, Caputo_Piecewise())
```

There many types of definitions of fractional derivative, Caputo is one of these useful definitions. The Caputo fractional derivative is first be proposed in [Michele Caputo's Paper](https://doi.org/10.1111/j.1365-246X.1967.tb02303.x), 

```math
^CD_t^\alpha f(t) = \frac{1}{\Gamma(n-\alpha)}\int_0^t\frac{f^{(n)}(\tau)d\tau}{(t-\tau)^{\alpha+1-n}}, n=\lceil{\alpha}\rceil
```

## Grünwald Letnikov sense derivative

```math
D^\alpha f(t)=\displaystyle \lim_{h\rightarrow0}\frac{1}{h^\alpha}\sum_{0\leq m\lt\infty}(-1)^m {{\alpha}\choose{m}}f(t+(\alpha-m)h)
```

## Highlight on some algorithms

FractionalCalculus.jl support many algorithms to calculate fractional derivative, here, I want to highlight the **Triangular Strip Matrix** method proposed by [Prof Igor](http://people.tuke.sk/igor.podlubny/index.html) to discrete fractional derivative.

```julia-repl
julia> fracdiff(f, α, end_point, h, RLDiff_Matrix())
```

With this advancing algorithms, we can not only compute the fractional derivative, but also the integer derivative! Or more precisely, arbitrary order!!!!🙌 Higher order derivative is also a piece of cake!!!!!!🎉

Try to set $\alpha$ as an integer, arbitrary integer of course! I promise you would enjoy it😏

