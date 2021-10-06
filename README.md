# FractionalCalculus.jl

<p align="center">
<img width="250px" src="https://raw.githubusercontent.com/ErikQQY/FractionalCalculus.jl/master/docs/src/assets/logo.png"/>
</p>

FractionalCalculus.jl provides support for fractional calculus computing.

## 🎇 Installation

If you have already install Julia, you can install FractionalCalculus.jl in REPL using Julia package manager:

```julia
import Pkg;
Pkg.add("FractionalCalculus")
```

Or if you want to experience the latest version of FractionalCalculus.jl, you can directly clone this repo:

```bash
git clone https://github.com/ErikQQY/FractionalCalculus.jl
```

## 🦸 Quick start

### Derivative

To compute the fractional derivative in a specific point, for example, compute <img src="https://latex.codecogs.com/gif.latex?\alpha=0.2" /> derivative of <img src="https://latex.codecogs.com/gif.latex?f(x)=x" /> in <img src="https://latex.codecogs.com/gif.latex?x=1" /> with step size <img src="https://latex.codecogs.com/gif.latex?h=0.0001" /> using **Caputo** sense:

```julia
Julia> fracdiff(x->x, 0.2, 0, 1, 0.0001, Caputo_Direct())
(1.0736712699679294, 1.4251042482525302e-8)
```

This will return a tuple **(result, estimating error)**.

### Integral

To compute the fractional integral in a specific point, for example, compute the semi integral of <img src="https://latex.codecogs.com/gif.latex?f(x)=x " /> in <img src="https://latex.codecogs.com/gif.latex?x=1" />  with step size <img src="https://latex.codecogs.com/gif.latex?h=0.0001" /> using **Riemann-Liouville** sense:

```julia
julia> fracint(x->x, 0.5, 0, 1, 0.0001, RL_Direct())
(0.7522527785271369, 8.022170098417246e-9)
```

This will return a tuple **(result, estimating error)**.

## 💻 All algorithms

```
Current Algorithms
├── FracDiffAlg
│   ├── Caputo
|   |   ├── Caputo_Direct
|   |   ├── Caputo_Direct_First_Diff_Known
|   |   ├── Caputo_Direct_First_Second_Diff_Known
|   |   └── Caputo_Piecewise
|   |
│   ├── GL
|   |   └── GL_Direct
|   |
|   └── RLDiff
|       └── RLDiff_Approx
|
└── FracIntAlg
    └── RLInt
        ├── RL_Direct
        ├── RL_Direct_First_Diff_Known
        ├── RL_Piecewise
        └── RLInt_Approx
```

## 📢 Status

Right now, FractionalCalculus.jl has only supports for little algorithms:

Fractional Derivative:

- [x] Caputo fractional derivative
- [x] Grunwald-Letnikov fractional derivative
- [ ] Caputo-Fabrizio fractional derivative
- [x] Riemann-Liouville fractional derivative 
- [ ] Atangana-Baleanu fractional derivative
- [ ] Riesz fractional derivative
- [ ] Marchaud fractional derivative
- [ ] Hadamard  fractional derivative
- [ ] Weyl  fractional derivative
- [ ] ......

Fractional Integral:
- [x] Riemann-Liouville fractional integral
- [ ] Hadamard fractional integral
- [ ] Atangana-Baleanu fractional integral
- [ ] ......

## 🧙 About Symbolic differentiation and integration

I am trying to find a way to support symbolic differentiation and integration features🤔.

## 🥂 Contributing

If you are interested in Fractional Calculus and Julia, welcome to raise an issue or file a Pull Request!!

This package is inspired by [JuliaMath/Calculus.jl](https://github.com/JuliMath/Calculus.jl) and [Numerical Methods in Fractional Calculus by Sean Townsend 2015](http://broncoscholar.library.cpp.edu/bitstream/10211.3/160926/1/TownsendSean_Thesis2015.pdf)