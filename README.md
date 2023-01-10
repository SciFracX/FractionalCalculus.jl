# FractionalCalculus.jl

<p align="center">
<img width="250px" src="https://raw.githubusercontent.com/SciFracX/FractionalCalculus.jl/master/docs/src/assets/logo.svg"/>
</p>


<p align="center">
  <a href="https://github.com/SciFracX/FractionalCalculus.jl/actions?query=workflow%3ACI">
    <img alt="building" src="https://github.com/SciFracX/FractionalCalculus.jl/workflows/CI/badge.svg">
  </a>
  <a href="https://codecov.io/gh/SciFracX/FractionalCalculus.jl">
    <img alt="codecov" src="https://codecov.io/gh/SciFracX/FractionalCalculus.jl/branch/master/graph/badge.svg">
  </a>
  <a href="https://scifracx.github.io/FractionalCalculus.jl/dev/">
    <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="license">
  </a>
  <a href="https://github.com/SciFracX/FractionalCalculus.jl/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/SciFracX/FractionalCalculus.jl?style=flat-square" alt="license">
  </a>
    <a href="https://zenodo.org/badge/latestdoi/420992306">
  	<img src="https://zenodo.org/badge/420992306.svg" alt="DOI">
  </a>
</p>

<p align="center">
  <a href="https://github.com/SciFracX/FractionalCalculus.jl/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/SciFracX/FractionalCalculus.jl?style=flat-square">
  </a>
  <a href="#">
    <img alt="GitHub stars" src="https://img.shields.io/github/stars/SciFracX/FractionalCalculus.jl?style=flat-square">
  </a>
  <a href="https://github.com/SciFracX/FractionalCalculus.jl/network">
    <img alt="GitHub forks" src="https://img.shields.io/github/forks/SciFracX/FractionalCalculus.jl?style=flat-square">
  </a>
</p>

FractionalCalculus.jl provides support for fractional calculus computing.

## 🎇 Installation

If you have already install Julia, you can install FractionalCalculus.jl in REPL using Julia package manager:

```julia
pkg> add FractionalCalculus
```

## 🦸 Quick start

### Derivative

To compute the fractional derivative in a specific point, for example, compute $ \alpha=0.2 $ derivative of $ f(x)=x $ in $ x=1 $ with step size $ h=0.0001 $ using **Riemann Liouville** sense:

```julia
julia> fracdiff(x->x, 0.2, 1, 0.0001, RLDiffL1())
1.0736712740308347
```

This will return the estimated value with high precision.

### Integral

To compute the fractional integral in a specific point, for example, compute the semi integral of $ f(x)=x $ in $ x=1 $  with step size $ h=0.0001 $ using **Riemann-Liouville** sense:

```julia
julia> fracint(x->x, 0.5, 1, 0.0001, RLIntApprox())
0.7522525439593486
```

This will return the estimated value with high precision.

## 💻 All algorithms

```
Current Algorithms
├── FracDiffAlg
│   ├── Caputo
|   |   ├── CaputoDirect
|   |   ├── CaputoPiecewise
|   |   ├── CaputoDiethelm
|   |   ├── CaputoHighPrecision
|   |   └── CaputoHighOrder
|   |
│   ├── Grünwald Letnikov
|   |   ├── GLDirect
|   |   ├── GLMultiplicativeAdditive
|   |   ├── GLLagrangeThreePointInterp
|   |   └── GLHighPrecision
|   |
|   ├── Riemann Liouville
|   |    ├── RLDiffL1
|   |    ├── RLLinearSplineInterp
|   |    ├── RLDiffMatrix
|   |    ├── RLG1
|   |    └── RLD
|   | 
|   ├── Hadamard
|   |    ├── HadamardLRect
|   |    ├── HadamardRRect
|   |    └── HadamardTrap
|   |
|   ├── Riesz
|   |    ├── RieszSymmetric
|   |    └── RieszOrtigueira
|   |
|   ├── Caputo-Fabrizio
|   |    └── CaputoFabrizioAS
|   |
|   └── Atanagana Baleanu
|        └── AtanganaSeda
|
└── FracIntAlg
    ├── Riemann Liouville
    |   ├── RLDirect
    |   ├── RLPiecewise
    |   ├── RLLinearInterp
    |   ├── RLIntApprox
    |   ├── RLIntMatrix
    |   ├── RLIntSimpson
    |   ├── RLIntTrapezoidal
    |   ├── RLIntRectangular
    |   └── RLIntCubicSplineInterp
    |
    └── Hadamard
        └── HadamardMat
```

For detailed usage, please refer to [our manual](https://scifracx.org/FractionalCalculus.jl/dev/Derivative/derivativeapi/).

## 🖼️ Example

Let's see examples here:

Compute the semi-derivative of $ f(x)=x $ in the interval $ [0, 1] $:

![Plot](/docs/src/assets/semiderivativeplot.png)

We can see the computing retains high precision⬆️.

Compute different order derivative of $ f(x)=x $:

![Different Order](/docs/src/assets/different_order_x_derivative.png)

Also different order derivative of $ f(x)=\sin(x) $:

![Different Order of sin](/docs/src/assets/different_order_sin_derivative.png)

And also different order integral of $ f(x)=x $:

![Different Order Of x](/docs/src/assets/different_order_x_integral.png)

<!---

Or arbitrary order derivative? A piece of cake!!😉

![Arbitrary](/docs/src/assets/arbitrary_order_derivative.png)

-->

## 🧙 Symbolic Fractional Differentiation and Integration

Thanks to SymbolicUtils.jl, FractionalCalculus.jl can do symbolic fractional differentiation and integration now!!

```julia
julia> using FractionalCalculus, SymbolicUtils
julia> @syms x
julia> semidiff(log(x))
log(4x) / sqrt(πx)
julia> semiint(x^4)
0.45851597901024005(x^4.5)
```

## 📢 Status

Right now, FractionalCalculus.jl has only supports for little algorithms:

Fractional Derivative:

- [x] Caputo fractional derivative
- [x] Grunwald-Letnikov fractional derivative
- [x] Riemann-Liouville fractional derivative 
- [x] Riesz fractional derivative
- [x] Hadamard  fractional derivative
- [x] Caputo-Fabrizio fractional derivative
- [x] Atangana-Baleanu fractional derivative
- [ ] Marchaud fractional derivative
- [ ] Weyl fractional derivative
- [ ] ......

Fractional Integral:
- [x] Riemann-Liouville fractional integral
- [x] Hadamard fractional integral
- [ ] Atangana-Baleanu fractional integral
- [ ] ......

## 📚 Reference

FractionalCalculus.jl is built upon the hard work of many scientific researchers, I sincerely appreciate what they have done to help the development of science and technology.

## 🥂 Contributing

If you are interested in Fractional Calculus and Julia, welcome to raise an issue or file a Pull Request!!

