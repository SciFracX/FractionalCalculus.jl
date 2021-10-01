module FractionalCalculus

using QuadGK
using SpecialFunctions: gamma

include("./Derivative/Derivative.jl")
include("./Integral/Integral.jl")

export FracDiffAlg, FracIntAlg
export fracdiff, fracint

#Export fractional derivative releating API
export Caputo, GL
export Caputo_Direct, Caputo_Direct_First_Diff_Known, Caputo_Direct_First_Second_Diff_Known

#Export fractional integral releating API
export RL
export RL_Direct, RL_Direct_First_Diff_Known, RL_Piecewise
end