# Fractional derivative

```@contents
Pages = ["derivative.md"]
```

To get start with fractional derivative, you need to know that unlike Newtonian derivatives, fractional derivative is defined via integral.

!!! tip "Non-Local Operators"
	It is noteworthy that the fractional derivatives are not local operators, which means that we cannot calculate the fractional derivative solely on the basis of function values of $f(x)$ taken from neighborhood of the point $x$. Instead, we have to take the entire history of $f(x)$ (i.e., all function values $f(x)$ for $0<x<a$) into account.

## Riemann Liouville sense derivative

Riemann Liouville sense derivative is built upon the Riemann Liouville sense integral.

```math
_aD^\alpha_tf(t)=\frac{d^n}{dt^n}\ _aD^{-(n-\alpha)}_tf(t)=\frac{d^n}{dt^n}\ _aI^{n-\alpha}_tf(t)
```

```math
_tD^\alpha_bf(t)=\frac{d^n}{dt^n}\ _tD^{-(n-\alpha)}_bf(t)=\frac{d^n}{dt^n}\ _tI^{n-\alpha}_bf(t)
```

We can use **FractionalCalculus.jl** to compute Riemann Liouville sense $0.5$ order fractional derivative of $f(x)=x$ at $x=1$ with step size $0.0001$:

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.0001, RLDiffL1())
1.1283791670955168
```

## Caputo sense derivative

There many types of definitions of fractional derivative, Caputo is one of these useful definitions. The Caputo fractional derivative is first be proposed in [Michele Caputo's Paper](https://doi.org/10.1111/j.1365-246X.1967.tb02303.x), 

```math
^CD_t^\alpha f(t) = \frac{1}{\Gamma(n-\alpha)}\int_0^t\frac{f^{(n)}(\tau)d\tau}{(t-\tau)^{\alpha+1-n}}, n=\lceil{\alpha}\rceil
```

In **FractionalCalculus.jl**, let's see, if you want to calculate the $0.5$ order fractional derivative of $f(x)=x$ at a $x=1$ with step size $0.0001$, simply typing these:


```julia-repl
julia>fracdiff(x->x, 0.5, 1, 0.0001, CaputoPiecewise())
1.128379167055761
```

We can see the result is closely resembling with the result in Riemann Liouville sense.

!!! info
	The Riemann Liouville sense derivative and Caputo sense derivative have the following relationship:
	```math
	_{RL}D^{\alpha}_{a, t}f(t)={_CD^\alpha_{a, t}f(t)}+\sum_{k=0}^{m-1}\frac{f^{(k)}(a)(t-a)^{k-a}}{\Gamma(k+1-\alpha)}
	```
	> Here ``m-1<\alpha<m``, ``f^{(m)}`` is integrable on ``[a, t]``.

## Grünwald Letnikov sense derivative

```math
D^\alpha f(t)=\displaystyle \lim_{h\rightarrow0}\frac{1}{h^\alpha}\sum_{0\leq m\lt\infty}(-1)^m {{\alpha}\choose{m}}f(t+(\alpha-m)h)
```

To compute the Grünwald Letnikov sense derivative, you can use FractionalCalculus.jl by:

```julia-repl
julia> fracdiff(x->x, 0.5, collect(0:0.01:1), 2, GLHighPrecision())
101-element Vector{Float64}:
 0.0
 0.11283791670955126
 0.15957691216057307
 0.19544100476116796
 0.22567583341910252
 0.252313252202016
 0.27639531957706837
 0.29854106607209235
 0.31915382432114614
 ⋮
 1.082303275961202
 1.0881694613449238
 1.094004191971426
 1.0998079684646789
 1.1055812783082735
 1.111324596323283
 1.1170383851240115
 1.1227230955528664
 1.1283791670955126
```

Here, we use the high precision algorithm, the fourth parameter means we set the precision order as **p=2**. The returned result means the derivative on the interval $[0, 1]$.

!!! info
	If the function ``f(t)`` is suitably smooth, then the Grünwald Letnikov sense derivative and the Riemann Liouville sense derivative is equivalent.

## Riesz sense derivative

The Riesz sense symmetric fractional derivative is defined by **Caputo** derivative:

```math

\frac{d^\beta \phi(x)}{d|x|^\mu}=D_{R}^{\beta}\phi(x)=\frac{1}{2}\Bigl({_{a}D}_{x}^{\beta}\phi(x)+{_{x}D_{b}^{\beta}}\phi(x)
 \Bigr)

```

In FractionalCalculus.jl, we can use the **RieszSymmetric** algorithm to compute the fractional derivative:

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.01, RieszSymmetric())
```

## Hadamard sense derivative

The Hadamard sense derivative is defined using Hadamard sense integral:

```math

{_HD_{a+}^{-\alpha}}f(x)=\frac{1}{\Gamma(\alpha)}\int_a^x(\log\frac{x}{t})^{-(1-\alpha)}f(t)\frac{dt}{t},\ x > a

```

So we can know the Hadamard sense fractional derivative:

```math

{_HD_{a+}^\alpha}f(x)=\delta^n[{_HD_{a+}^{-(n-\alpha)}f(x)}],\ x>a

```

```math

\delta=x\frac{d}{dx},\ n-1\leq\alpha<n\in\mathrm{Z^+}

```

To compute Hadamard fractional derivative, we can use the Hadamard relating algorithms in FractionalCalculus.jl:

```julia-repl
julia> fracdiff(x->x, 0.5, 0, 1, 0.01, HadamardLRect())
0.9165222777761635
```

!!! note "Non-singular kernel"
	When we are using the "fractional derivative" with non-singular kernel, like Caputo-Fabrizio and Atangana-Baleanu-Caputo sense, we should know these operators are in fact integer order operators, for more details we recommend you to read the paper [Why Fractional Derivatives with Nonsingular Kernels Should Not Be Used](https://link.springer.com/article/10.1515/fca-2020-0032)

## Caputo-Fabrizio sense derivative

Caputo-Fabrizio sense fractional order derivative is defined by:

```math
{^{CF}_0D^\alpha_t u(t)}=\frac{M(\alpha)}{1-\alpha}\int^t_0\frac{d}{d\tau}u(\tau)\exp[-\frac{\alpha}{1-\alpha}(t-\tau)]d\tau\\

M(\alpha)=1-\alpha+\alpha/\Gamma(\alpha)
```

To compute the Caputo-Fabrizio sense derivative, we can use the ```CaputoFabrizioAS``` algorithm in FractionalCalculus.jl:

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.01, CaputoFabrizioAS())
0.9887564512257243
```

## Atangana-Baleanu sense derivative

Atangana-Baleanu sense fractional order derivative is defined by:

```math
{^{ABC}_0D^\alpha_t f(t)}=\frac{AB(\alpha)}{1-\alpha}\int^t_0\frac{d}{d\tau}f(\tau)E_\alpha[-\frac{\alpha}{1-\alpha}(t-\tau)^\alpha]d\tau\\

AB(\alpha)=1-\alpha+\alpha/\Gamma(\alpha)
```

To compute the Atangana-Baleanu sense derivative, we can use the ```AtanganaSeda``` algorithm in FractionalCalculus.jl:

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.01, AtanganaSeda())
-0.8696378200415389
```


!!! note
	Here we need to specify the **start point** and **end point**

There are different approximating methods being used in the computing, choose the one you need😉