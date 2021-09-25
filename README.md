# FractionalCalculus.jl

![Logo](/docs/logo.png)

This package provides support for fractional calculus computing.

## Installation

If you have already install Julia, you can install FractionalCalculus.jl in REPL using Julia package manager:

```julia
import Pkg;
Pkg.add("FractionalCalculus")
```

Or if you want to experience the latest version of FractionalCalculus.jl, you can directly clone this repo:

```bash
git clone https://github.com/ErikQQY/FractionalCalculus.jl
```

## Status

Right now, FractionalCalculus.jl has only supports for little algorithms:

Fractional Derivative:

- [x] Caputo fractional derivative
- [x] Grunwald-Letnikov fractional derivative
- [ ] Caputo-Fabrizio fractional derivative
- [ ] Riemann-Liouville fractional derivative 
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

## About Symbolic differentiation and integration

Trying to find a way to support symbolic differentiation and integration features

## Contributing

If you are interested in Fractional Calculus and Julia, welcome to raise an issue or file a Pull Request!!

This package is inspired by [JuliaMath/Calculus.jl](https://github.com/JuliMath/Calculus.jl) and [Numerical Methods in Fractional Calculus by Sean Townsend 2015](http://broncoscholar.library.cpp.edu/bitstream/10211.3/160926/1/TownsendSean_Thesis2015.pdf)