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
```

```@docs
FractionalCalculus.Caputo_Direct
FractionalCalculus.Caputo_Direct_First_Diff_Known
FractionalCalculus.Caputo_Direct_F_Second_Diff_Known
FractionalCalculus.Caputo_Piecewise
FractionalCalculus.Caputo_Diethelm
FractionalCalculus.GL_Direct
FractionalCalculus.GL_Finite_Difference
FractionalCalculus.GL_Multiplicative_Additive
FractionalCalculus.GL_Lagrange_Three_Point_Interp
FractionalCalculus.RLDiff_Approx
FractionalCalculus.RLDiff_Matrix
```

```@docs
FractionalCalculus.@fracdiff
FractionalCalculus.@semifracdiff
```

