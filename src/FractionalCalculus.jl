module FractionalCalculus

using QuadGK
using SpecialFunctions: gamma
using LinearAlgebra, InvertedIndices
using SpecialMatrices: Vandermonde
using SymbolicUtils, SymbolicUtils.Rewriters
using Symbolics

include("./Derivative/Derivative.jl")
include("./Integral/Integral.jl")
include("./Derivative/mlfun.jl")

export FracDiffAlg, FracIntAlg
export fracdiff, fracint

# Export symbolics computing APIs
export semidiff, semiint
export SEMIDIFFRULES, SEMIINTRULES
export Incomplete_beta, AuxiliaryFresnelSin, AuxiliaryFresnelCos, Struve, MStruve, Legendre, SinIntegral, HyperSinIntegral, Hypergeometric1F1

# Export fractional derivative releating API
export Caputo, GL, RLDiff, Hadamard

# Caputo sense fractional derivative
export CaputoDirect, CaputoPiecewise, CaputoDiethelm, CaputoHighPrecision, CaputoHighOrder

# Grünwald Letnikov sense fractional derivative
export GLDirect, GLMultiplicativeAdditive, GLLagrangeThreePointInterp, GLFiniteDifference, GLHighPrecision

# Riemann Liouville sense fractional derivative
export RLDiffApprox, RLDiffMatrix, RLLinearSplineInterp, RLG1, RLD

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
export RieszMatrix, omega, B, genfun, first_order

export mittleff

end
