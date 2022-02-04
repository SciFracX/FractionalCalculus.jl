import FractionalCalculus: FracDiffAlg, omega, B

"""
Riesz sense fractional derivative
"""
abstract type Riesz <: FracDiffAlg end

"""
# Riesz sense symmetric fractional derivative algorithm.

    fracdiff(f, α, end_point, h, Riesz_Symmetric())

Compute fractional derivative of Riesz sense using Triangular Strip Matrix algorithm.

### Example

```julia-repl
julia> fracdiff(x->x, 0.5, 1, 0.01, Riesz_Symmetric())
```
"""
struct Riesz_Symmetric <: Riesz end


################################################################
###                    Type definition done                  ###
################################################################


function fracdiff(f, α, end_point, h, ::Riesz_Symmetric)
    N=Int(floor(end_point/h))

    mat = RieszMatrix(α, N+1, h)
    return mat*f.(collect(0:h:end_point))
end
function RieszMatrix(α, N, h)
    caputo = B(N+1, α)
    caputo = caputo[2:(N+1), 1:N]
    result = 1/2*(caputo+caputo')
    result = h^(-α)*result

    return result
end