"""
"""
abstract type ChebyshevDiff <: FracDiffAlg end

@doc raw"""

Compute ``\partial^{\alpha} f(s) `` for ``0 < \alpha < 1``.
"""
struct ChebyshevSeriesUnit <: ChebyshevDiff end

function chebyshev_U(t::K, n::Integer) where K <: Number
    U::K = if n == 0
        1
    elseif n == 1
        2t
    else
        Uᵢ = 1
        Uⱼ = 2t

        for _ = 1:n-1
            Uᵢ, Uⱼ = (Uⱼ, 2t * Uⱼ - Uᵢ)
        end

        Uⱼ
    end

    return U
end

function fracdiff(f::Union{Function, K}, α::K, s::K, n::Integer, ::ChebyshevSeriesUnit) where K <: Number
    @assert 0 < α < 1

    # ~ Aliasing ~ #
    U = chebyshev_U

    # ~ Step 1: f(x); given
    F = f isa K ? (::K) -> f : f

    # ~ Step 2: Compute coefficients of expansion into the Chebyshev series; {aₖ}
    a = zeros(n + 1)

    for k = 0:n
        τₖ::K = if 1 < k < n
            2k / n
        else
            1 / 2n
        end
        
        a[k + 1] = τₖ * sum( F(2 * cos(π * l / n) - 1) * cos(π * k * l / n) for l = 0:n )
    end

    # ~ Step 3: Compute coefficients of expansion into the Chebyshev series; {bₖ}
    b = zeros(K, n + 1)

    for k = n-1:1
        b[k] = ((4s - 2) * b[k + 1] - (1 - (1 - α) / (k + 2)) * b[k + 2] + (8k + 8) * a[k + 2]) / (1 + (1 - α) / k)
    end

    # ~ Step 4: Compute the value of fₙ'(s)
    f_s::K = 2 * sum((k + 1) * a[k + 2] * U(2s - 1, k + 1) for k = 0:n-1)

    # ~ Step 5: Compute the value of F(x) at x = 0 and x = s
    c = K[((b[k] / 4k) - (b[k + 2] / (4k + 2))) for k = 1:n-1]

    F_0::K = sum(c[k] * U(-1, k) for k = 1:n-1)
    F_s::K = sum(c[k] * U(2s - 1, k) for k = 1:n-1)

    # ~ Step 6: Compute the value of ∂ᵅf(s)
    return F(0) * s ^ (-α) + (f_s / (1 - α) - F_s + F_0) * s ^ (1 - α)
end


@doc raw"""
"""
function chebyshev_trait(α::T) where T <: Number
    if 0 < α < 1
        return ChebyshevSeriesUnit()
    else
        error("No algorithms for α = $(α)")
    end
end

struct ChebyshevSeries <: ChebyshevDiff end

fracdiff(f, α, s, n, ::ChebyshevSeries) = fracdiff(f, α, s, n, chebyshev_trait(α))