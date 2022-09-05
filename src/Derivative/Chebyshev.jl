"""
"""
abstract type ChebyshevDiff <: FracDiffAlg end

@doc raw"""

Compute ``\partial^{\alpha} f(s) `` for ``0 < \alpha < 1``.
"""
struct ChebyshevSeriesUnit <: ChebyshevDiff end

function _chebyshev_U(t::K, n::Integer) where {K<:Number}
    if n == 1     # k = 0
        K[1]
    elseif n == 2 # k = 0, 1
        K[1, 2t]
    else          # k = 0, 1, 2 ... n - 1
        U = Vector{K}(undef, n)

        U[1] = 1
        U[2] = 2t

        for k = 3:n
            U[k] = 2t * U[k-1] - U[k-2]
        end

        return U
    end
end

function fracdiff(f::Union{Function,K}, α::K, s::K, n::Integer, ::ChebyshevSeriesUnit) where {K<:Number}
    @assert 0 < α < 1

    # ~ Aliasing ~ #
    U = _chebyshev_U

    # ~ Step 1: f(x); given
    F = f isa K ? (::K) -> f : f

    # ~ Step 2: Compute coefficients of expansion into the Chebyshev series; {aₖ}
    a = zeros(K, n + 1)

    for k = 0:n
        τ = 0 < k < n ? 2k : 1
        ρ = 0

        for l = 0:n
            ρ += if l == 0
                F(1) / 2
            elseif l == n
                (F(-3) * (-1)^k) / 2
            else
                F(2 * cos(π * l / n) - 1) * cos((π * k * l) / n)
            end
        end

        a[k+1] = (τ / n) * ρ
    end

    # ~ Step 3: Compute coefficients of expansion into the Chebyshev series; {bₖ}
    b = zeros(K, n + 1)

    for k = (n-1):-1:1
        b[k] = (
            2(2s - 1) * b[k+1]
            - (1 - (1 - α) / (k + 2)) * b[k+2]
            + 8(k + 1) * a[k+2]
        ) / (1 + (1 - α) / k)
    end

    # ~ Compute Chebyshev Polynomials
    U_0 = U(   - 1, n) # ~ x = 0; k = 0, 1, ..., n - 1
    U_s = U(2s - 1, n) # ~ x = s; k = 0, 1, ..., n - 1

    # ~ Step 4: Compute the value of fₙ'(s)
    df_s = 2 * sum((k + 1) * a[k+2] * U_s[k+1] for k = 0:n-1)

    # ~ Step 5: Compute the value of F(x) at x = 0 and x = s
    c = [(b[k] / 4k) - (b[k+2] / 4(k + 2)) for k = 1:n-1]

    F_0 = sum([c[k] * U_0[k+1] for k = 1:n-1])
    F_s = sum([c[k] * U_s[k+1] for k = 1:n-1])

    # ~ Step 6: Compute the value of ∂ᵅf(s)
    return F(0) * s^(-α) + ((df_s / (1 - α)) - F_s + F_0) * s^(1 - α)
end


@doc raw"""
"""
function chebyshev_trait(α::T) where {T<:Number}
    if 0 < α < 1
        return ChebyshevSeriesUnit()
    else
        error("No algorithms for α = $(α)")
    end
end

struct ChebyshevSeries <: ChebyshevDiff end

fracdiff(f, α, s, n, ::ChebyshevSeries) = fracdiff(f, α, s, n, chebyshev_trait(α))