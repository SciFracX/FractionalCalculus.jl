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

> See our talk on JuliaCN 2021 Winter Conf: [Slide](https://julia-cn-conf2021.vercel.app/1), [Youtube](https://www.youtube.com/watch?v=oVvrW7EgEwg), [BiliBili](https://www.bilibili.com/video/BV1vY411W7Dw?p=18)

FractionalCalculus.jl provides support for fractional calculus computing.

## 🎇 Installation

If you have already install Julia, you can install FractionalCalculus.jl in REPL using Julia package manager:

```julia
Pkg> add FractionalCalculus
```

Or if you want to experience the latest version of FractionalCalculus.jl:

```julia
Pkg> add FractionalCalculus#master
```

## 🦸 Quick start

### Derivative

To compute the fractional derivative in a specific point, for example, compute <img src="https://latex.codecogs.com/svg.image?\inline&space;\alpha=0.2" title="\inline \alpha=0.2" /> derivative of <img src="https://latex.codecogs.com/svg.image?\inline&space;f(x)=x" title="\inline f(x)=x" /> in <img src="https://latex.codecogs.com/svg.image?\inline&space;x=1" title="\inline x=1" /> with step size <img src="https://latex.codecogs.com/svg.image?\inline&space;h=0.0001" title="\inline h=0.0001" /> using **Riemann Liouville** sense:

```julia
julia> fracdiff(x->x, 0.2, 1, 0.0001, RLDiff_Approx())
1.0736712740308347
```

This will return the estimated value with high precision.

### Integral

To compute the fractional integral in a specific point, for example, compute the semi integral of <img src="https://latex.codecogs.com/svg.image?\inline&space;f(x)=x" title="\inline f(x)=x" /> in <img src="https://latex.codecogs.com/svg.image?\inline&space;x=1" title="\inline x=1" />  with step size <img src="https://latex.codecogs.com/svg.image?\inline&space;h=0.0001" title="\inline h=0.0001" /> using **Riemann-Liouville** sense:

```julia
julia> fracint(x->x, 0.5, 1, 0.0001, RLInt_Approx())
0.7522525439593486
```

This will return the estimated value with high precision.

## 💻 All algorithms

```
Current Algorithms
├── FracDiffAlg
│   ├── Caputo
|   |   ├── Caputo_Direct
|   |   ├── Caputo_Direct_First_Diff_Known
|   |   ├── Caputo_Direct_First_Second_Diff_Known
|   |   ├── Caputo_Piecewise
|   |   ├── Caputo_Diethelm
|   |   └── Caputo_High_Precision
|   |
│   ├── Grünwald Letnikov
|   |   ├── GL_Direct
|   |   ├── GL_Multiplicative_Additive
|   |   ├── GL_Lagrange_Three_Point_Interp
|   |   └── GL_High_Precision
|   |
|   ├── Riemann Liouville
|   |    ├── RLDiff_Approx
|   |    ├── RL_Linear_Spline_Interp
|   |    └── RLDiff_Matrix
|   | 
|   ├── Hadamard
|   |    ├── Hadamard_LRect
|   |    ├── Hadamard_RRect
|   |    └── Hadamard_Trap
|   |
|   └── Riesz
|        └── Riesz_Symmetric
|
└── FracIntAlg
    └── Riemann Liouville
        ├── RL_Direct
        ├── RL_Direct_First_Diff_Known
        ├── RL_Piecewise
        ├── RL_LinearInterp
        ├── RLInt_Approx
        ├── RLInt_Matrix
        ├── RLInt_Simpson
        ├── RLInt_Trapezoidal
        ├── RLInt_Rectangular
        └── RLInt_Cubic_Spline_Interp
```

## 🖼️ Example

Let's see examples here:

Compute the semi-derivative of <img src="https://latex.codecogs.com/svg.image?\inline&space;f(x)=x" title="\inline f(x)=x" /> in the interval [0, 1]:

![Plot](/docs/src/assets/semiderivativeplot.png)

We can see the computing retains high precision⬆️.

Compute different order derivative of <img src="https://latex.codecogs.com/svg.image?\inline&space;f(x)=x" title="\inline f(x)=x" />:

![Different Order](/docs/src/assets/different_order_x_derivative.png)

Also different order derivative of <img src="https://latex.codecogs.com/svg.image?\inline&space;f(x)=\sin(x)" title="\inline f(x)=\sin(x)" />:

![Different Order of sin](/docs/src/assets/different_order_sin_derivative.png)

And also different order integral of <img src="https://latex.codecogs.com/svg.image?\inline&space;f(x)=x" title="\inline f(x)=x" />:

![Different Order Of x](/docs/src/assets/different_order_x_integral.png)

Or arbitrary order derivative? A piece of cake!!😉

![Arbitrary](/docs/src/assets/arbitrary_order_derivative.png)

## 📢 Status

Right now, FractionalCalculus.jl has only supports for little algorithms:

Fractional Derivative:

- [x] Caputo fractional derivative
- [x] Grunwald-Letnikov fractional derivative
- [x] Riemann-Liouville fractional derivative 
- [x] Riesz fractional derivative
- [x] Hadamard  fractional derivative
- [ ] Caputo-Fabrizio fractional derivative
- [ ] Atangana-Baleanu fractional derivative
- [ ] Marchaud fractional derivative
- [ ] Weyl  fractional derivative
- [ ] ......

Fractional Integral:
- [x] Riemann-Liouville fractional integral
- [ ] Hadamard fractional integral
- [ ] Atangana-Baleanu fractional integral
- [ ] ......

## 🧙 About Symbolic differentiation and integration

I am trying to find a way to support symbolic differentiation and integration features🤔.

## 📚 Reference

FractionalCalculus.jl is built upon the hard work of many scientific researchers, I sincerely appreciate what they have done to help the development of science and technology.

## 🥂 Contributing

If you are interested in Fractional Calculus and Julia, welcome to raise an issue or file a Pull Request!!

