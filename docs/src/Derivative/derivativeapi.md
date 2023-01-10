# Fractional derivative API

FractionalCalculus.jl has a disciplinary type system, by specifying which type you want to use, you can use the relating algorithm to compute fractional differentiation and fractional integral. 

This page contains all of the existing API we can use for computing fractional derivative.

```@docs
FractionalCalculus.fracdiff
```

## Caputo Sense fractional derivative

```@docs
FractionalCalculus.CaputoDirect
FractionalCalculus.CaputoPiecewise
FractionalCalculus.CaputoDiethelm
FractionalCalculus.CaputoHighPrecision
FractionalCalculus.CaputoHighOrder
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
FractionalCalculus.L1
FractionalCalculus.RLDiffMatrix
FractionalCalculus.RLG1
FractionalCalculus.RLD
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