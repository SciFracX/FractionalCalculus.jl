module FractionalCalculus

using QuadGK
using SpecialFunctions: gamma
using LinearAlgebra, InvertedIndices, SpecialMatrices

include("./Derivative/Derivative.jl")
include("./Integral/Integral.jl")

export FracDiffAlg, FracIntAlg
export fracdiff, fracint

# Export fractional derivative releating API
export Caputo, GL, RLDiff, Hadamard

# Caputo sense fractional derivative
export Caputo_Direct, Caputo_Direct_First_Diff_Known, Caputo_Direct_First_Second_Diff_Known
export Caputo_Piecewise
export Caputo_Diethelm
export GL_Direct

# Grünwald Letnikov sense fractional derivative
export GL
export GL_Multiplicative_Additive, GL_Lagrange_Three_Point_Interp

# Riemann Liouville sense fractional derivative
export RLDiff_Approx
export RLDiff_Matrix

# Hadamard sense fractional derivative
export Hadamard_LRect, Hadamard_RRect, Hadamard_Trap

# Riesz sense fractional derivative
export Riesz_Symmetric

#Export fractional integral relating API
export RLInt
export RL_Direct, RL_Direct_First_Diff_Known, RL_Piecewise
export RLInt_Approx, RL_LinearInterp, RLInt_Matrix
end