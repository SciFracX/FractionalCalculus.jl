import FractionalCalculus.FracDiffAlg

"""
Riesz sense fractional derivative algorithms.
"""
abstract type Riesz <: FracDiffAlg end

"""
Using the Triangular Strip Matrices to compute the fractional derivative Riesz derivative.
"""
struct Riesz_Symmetric <: Riesz end


################################################################
###                    Type defination done                  ###
################################################################


"""
    fracdiff(f, α, end_point, h, Riesz_Symmetric())

Compute fractional derivative of Riesz sense using Triangular Strip Matrix algorithm.
"""
function fracdiff(f, α, end_point, h, ::Riesz_Symmetric)
    N=Int(floot(end_point/h))

    mat = RieszMatrix(α, N+1, h)
    return mat*f.(collect(0:h:end_point))
end
function RieszMatrix(α, N, h)
    caputo=B(N+1, α)
    caputo=caputo[2:(N+1), 1:N]
    result=1/2*(caputo+caputo')
    result=h^(-α)*result

    return result
end
