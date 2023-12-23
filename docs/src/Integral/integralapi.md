# Fractional integral API

FractionalCalculus.jl has a disciplinary type system, by specifying which type you want to use, you can use the related algorithm to compute fractional differentiation and fractional integral. 

This page contains all of the existing API we can use for computing fractional integral.

```@docs
FractionalCalculus.fracint
```

## Riemann Liouville sense fractional integral

```@docs
FractionalCalculus.RLInt
FractionalCalculus.RLDirect
FractionalCalculus.RLPiecewise
FractionalCalculus.RLIntApprox
FractionalCalculus.RLLinearInterp
FractionalCalculus.RLIntMatrix
FractionalCalculus.RLIntSimpson
FractionalCalculus.RLIntTrapezoidal
FractionalCalculus.RLIntRectangular
FractionalCalculus.RLIntCubicSplineInterp
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