import FractionalCalculus.FracDiffAlg

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
"""
struct Riesz_Symmetric <: Riesz end


################################################################
###                    Type defination done                  ###
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

function B(N, p)
    result=zeros(N, N)
    temp=omega(N, p)

    @inbounds @simd for i ∈ 1:N
        @views result[i, 1:i]=reverse(temp[1:i])
    end

    return result
end

function omega(n, p)
    omega = zeros(n+1)

    omega[1]=1
    @fastmath @inbounds @simd for i ∈ 1:n
        omega[i+1]=(1-(p+1)/i)*omega[i]
    end
    
    return omega
end