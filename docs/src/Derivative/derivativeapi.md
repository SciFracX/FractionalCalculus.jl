# Fractional derivative API

FractionalCalculus.jl has a disciplinary type system, by specifying which type you want to use, you can use the related algorithm to compute fractional differentiation and fractional integral. 

This page contains all of the existing APIs we can use for computing fractional derivatives.

```@docs
FractionalCalculus.fracdiff
```

## Caputo Sense fractional derivative

```@docs
FractionalCalculus.CaputoDirect
FractionalCalculus.CaputoTrap
FractionalCalculus.CaputoDiethelm
FractionalCalculus.CaputoHighPrecision
FractionalCalculus.CaputoHighOrder
FractionalCalculus.CaputoL1
FractionalCalculus.CaputoL2
```

## Grunwald Letnikov sense fractional derivative

```@docs
FractionalCalculus.GLDirect
FractionalCalculus.GLFiniteDifference
FractionalCalculus.GLMultiplicativeAdditive
FractionalCalculus.GLLagrangeThreePointInterp
```

## Riemann Liouville sense fractional derivative

```@docs
FractionalCalculus.RLDiffL1
FractionalCalculus.RLDiffMatrix
FractionalCalculus.RLLinearSplineInterp
FractionalCalculus.RLG1
FractionalCalculus.RLD
FractionalCalculus.RLDiffL2
FractionalCalculus.RLDiffL2C
```

## Hadamard sense fractional derivative

```@docs
FractionalCalculus.HadamardLRect
FractionalCalculus.HadamardRRect
FractionalCalculus.HadamardTrap
```

## Riesz sense fractional derivative

```@docs
FractionalCalculus.RieszSymmetric
FractionalCalculus.RieszOrtigueira
```

## Atangana Baleanu sense fractional derivative

```@docs
FractionalCalculus.AtanganaSeda
```

## Caputo Fabrizio sense fractional derivative

```@docs
FractionalCalculus.CaputoFabrizioAS
```



```@docs
FractionalCalculus.@fracdiff
FractionalCalculus.@semifracdiff
```

## Symbolic fractional differentiation

```@docs
FractionalCalculus.semidiff
```