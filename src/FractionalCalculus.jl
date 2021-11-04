module FractionalCalculus

using QuadGK
using SpecialFunctions: gamma
using LinearAlgebra, InvertedIndices

include("./Derivative/Derivative.jl")
include("./Integral/Integral.jl")

export FracDiffAlg, FracIntAlg
export fracdiff, fracint

#Export fractional derivative releating API
export Caputo, GL, RLDiff
export Caputo_Direct, Caputo_Direct_First_Diff_Known, Caputo_Direct_First_Second_Diff_Known
export Caputo_Piecewise
export Caputo_Diethelm
export GL_Direct

export GL_Multiplicative_Additive, GL_Lagrange_Three_Point_Interp
export RLDiff_Approx

#Export fractional integral relating API
export RLInt
export RL_Direct, RL_Direct_First_Diff_Known, RL_Piecewise
export RLInt_Approx, RL_LinearInterp
end