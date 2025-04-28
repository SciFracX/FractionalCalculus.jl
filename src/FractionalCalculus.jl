module FractionalCalculus

using QuadGK
using SpecialFunctions
using SpecialFunctions: gamma
using LinearAlgebra, InvertedIndices
using SpecialMatrices: Vandermonde
using ForwardDiff
using UnPack

include("./Derivative/Derivative.jl")
include("./Integral/Integral.jl")
include("./Derivative/mlfun.jl")

export FracDiffAlg, FracIntAlg
export fracdiff, fracint

# Export symbolics computing APIs
"""
# Symbolic fractional differentiation

    semidiff(fun)

```semidiff``` uses SymbolicUtils.jl and Symbolics.jl to compute symbolic fractional differentiation.

### Example

```julia-repl
julia> using SymbolicUtils
julia> @syms x
julia> semidiff(log(x))
log(4x) / sqrt(πx)
```
"""
function semidiff end

"""
# Symbolic fractional integral

    semiint(fun)

```semiint``` uses SymbolicUtils.jl and Symbolics.jl to compute symbolic fractional integral.

### Example

```julia-repl
julia> using SymbolicUtils
julia> @syms x
julia> semiint(x^4)
0.45851597901024005(x^4.5)
```
"""
function semiint end

export semidiff, semiint

# Export fractional derivative releating API
export Caputo, GL, RLDiff

# Caputo sense fractional derivative
export CaputoDirect, CaputoTrap, CaputoDiethelm, CaputoHighPrecision, CaputoHighOrder, CaputoL1, CaputoL2

# Grünwald Letnikov sense fractional derivative
export GLDirect, GLMultiplicativeAdditive, GLLagrangeThreePointInterp, GLFiniteDifference, GLHighPrecision

# Riemann Liouville sense fractional derivative
export RLDiffL1, RLDiffL2, RLDiffL2C, RLDiffMatrix, RLLinearSplineInterp, RLG1, RLD

# Hadamard sense fractional derivative
export HadamardLRect, HadamardRRect, HadamardTrap

# Riesz sense fractional derivative
export RieszSymmetric, RieszOrtigueira

export AtanganaBaleanu, AtanganaSeda

# Caputo-Fabrizio sense fractional derivative
export CaputoFabrizio, CaputoFabrizioAS

# Macros for convenient computing
export @fracdiff, @semifracdiff

#Export fractional integral relating API
export RLInt, HadamardInt


export RLDirect, RLPiecewise
export RLIntApprox, RLLinearInterp, RLIntMatrix
export RLIntSimpson, RLIntTrapezoidal, RLIntRectangular, RLIntCubicSplineInterp

export HadamardMat

# Macros for convenient computing
export @fracint, @semifracint

# Auxiliary functions
export RieszMatrix, B, genfun

export mittleff

end
