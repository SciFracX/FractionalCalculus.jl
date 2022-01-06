# Fractional derivative API

FractionalCalculus.jl has a disciplinary type system, by specifying which type you want to use, you can use the relating algorithm to compute fractional differentiation and fractional integral. 

This page contain all of the existing API we can use for computing fractional derivative.

```@docs
FractionalCalculus.fracdiff
```

```@docs
FractionalCalculus.FracDiffAlg
```

```@docs
FractionalCalculus.Caputo
FractionalCalculus.GL
FractionalCalculus.RLDiff
FractionalCalculus.Hadamard
FractionalCalculus.Riesz
```

## Caputo Sense fractional derivative

```@docs
FractionalCalculus.Caputo_Direct
FractionalCalculus.Caputo_Piecewise
FractionalCalculus.Caputo_Diethelm
FractionalCalculus.Caputo_High_Precision
```

## Grunwald Letnikov sense fractional derivative

```@docs
FractionalCalculus.GL_Direct
FractionalCalculus.GL_Finite_Difference
FractionalCalculus.GL_Multiplicative_Additive
FractionalCalculus.GL_Lagrange_Three_Point_Interp
```

## Riemann Liouville sense fractional derivative

```@docs
FractionalCalculus.RLDiff_Approx
FractionalCalculus.RLDiff_Matrix
```

## Hadamard sense fractional derivative

```@docs
FractionalCalculus.Hadamard_LRect
FractionalCalculus.Hadamard_RRect
FractionalCalculus.Hadamard_Trap
```

## Riesz sense fractional derivative

```@docs
FractionalCalculus.Riesz_Symmetric
```

```@docs
FractionalCalculus.@fracdiff
FractionalCalculus.@semifracdiff
```

