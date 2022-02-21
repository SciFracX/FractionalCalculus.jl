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


"""
Ortigueira's definition of the symmetric Riesz derivative via centred differences
"""
struct Riesz_Ortigueira <: Riesz end

################################################################
###                    Type definition done                  ###
################################################################


function fracdiff(f::FunctionAndNumber, α, end_point, h, ::Riesz_Symmetric)
    N=floor(Int, end_point/h)

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

#FIXME: Values are changing when h become smaller?
function fracdiff(f, α, point, h, ::Riesz_Ortigueira)
    N = round(Int, point/h)
    t=collect(0:h:point)
    return ranort(α, N+1, h)*f.(t)
end

function ranort(alpha, N)
    k=collect(0:N-1)
    rc = ((-1)*ones(size(k))).^k.*gamma.(alpha+1).*(gamma.(alpha*0.5 .-k.+1).*gamma.(alpha*0.5 .+ k.+1)).^(-1)
    rc = rc*(cos(alpha*π*0.5))

    R = zeros(N, N)

    for m=1:N
        R[m, m:N] = rc[1:N-m+1]
    end

    for i=1:N-1
        for j=i:N
            R[j, i] = R[i, j]
        end
    end
    return R
end
# Multiple dispatch for ranort
function ranort(alpha, N, h)
    R = ranort(alpha, N)
    R = R*h^(-alpha)
    return R
end