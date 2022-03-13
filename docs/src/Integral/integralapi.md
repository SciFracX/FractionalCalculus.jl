# Fractional integral API

FractionalCalculus.jl has a disciplinary type system, by specifying which type you want to use, you can use the relating algorithm to compute fractional differentiation and fractional integral. 

This page contains all of the existing API we can use for computing fractional integral.

```@docs
FractionalCalculus.fracint
```

## Riemann Liouville sense fractional integral

```@docs
FractionalCalculus.RLInt
FractionalCalculus.RL_Direct
FractionalCalculus.RL_Piecewise
FractionalCalculus.RLInt_Approx
FractionalCalculus.RL_LinearInterp
FractionalCalculus.RLInt_Matrix
FractionalCalculus.RLInt_Simpson
FractionalCalculus.RLInt_Trapezoidal
FractionalCalculus.RLInt_Rectangular
FractionalCalculus.RLInt_Cubic_Spline_Interp
```

## Hadamard sense fractional integral

```@dcos
FractionalCalculus.HadamardMat
```

```@docs
FractionalCalculus.@fracint
FractionalCalculus.@semifracint
```

## Symbolic fractional integral

```@docs
FractionalCalculus.semiint
```