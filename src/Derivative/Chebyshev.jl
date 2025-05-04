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

function _chebyshev_a(::Type{K}, f::Function, n::Integer) where {K<:Number}
    a = Vector{K}(undef, n + 1)

    for k = 0:n
        τ = 1 < k < n ? 2k : 1
        ρ = zero(K)

        for l = 0:n
            x::K = 2 * cos(π * l / n) - 1
            y::K = cos((π * k * l) / n)
            z::K = f(x) * y

            ρ += 1 < l < n ? z : z / 2
        end

        a[k+1] = (τ / n) * ρ
    end

    return a
end

function _chebyshev_b(α::K, s::K, a::Vector{K}, n::Integer) where {K<:Number}
    @assert length(a) == n + 1

    b = Vector{K}(undef, n + 1)

    b[n+1] = b[n] = 0

    for k = (n-1):-1:1
        x::K = 2(2s - 1) * b[k+1]
        y::K = (1 - (1 - α) / (k + 2)) * b[k+2]
        z::K = 8(k + 1) * a[k+2]
        w::K = (1 + (1 - α) / k)

        b[k] = (x - y + z) / w
    end

    return b
end

function _chebyshev_c(b::Vector{K}, n::Integer) where {K<:Number}
    return K[(b[k] / 4k) - (b[k+2] / 4(k + 2)) for k = 1:n-1]
end

function fracdiff(_f::Union{Function,K}, α::K, s::K, n::Integer, ::ChebyshevSeriesUnit) where {K<:Number}
    @assert 0 < α < 1

    # ~ Step 1: f(x); given
    f = _f isa K ? (::K) -> _f : _f

    # ~ Step 2: Compute coefficients of expansion into the Chebyshev series; {aₖ}
    a = _chebyshev_a(K, f, n)

    # ~ Step 3: Compute coefficients of expansion into the Chebyshev series; {bₖ}
    b = _chebyshev_b(α, s, a, n)

    # ~ Compute Chebyshev Polynomials
    U_0 = _chebyshev_U(-1, n)     # Uₖ(2x - 1) ~ x = 0; k = 0, 1, ..., n - 1
    U_s = _chebyshev_U(2s - 1, n) # Uₖ(2x - 1) ~ x = s; k = 0, 1, ..., n - 1

    # ~ Step 4: Compute the value of fₙ'(s)
    df_s = 2 * sum(k * a[k+1] * U_s[k] for k = 1:n)

    # ~ Step 5: Compute the value of F(x) at x = 0 and x = s
    #      bₖ   bₖ₊₂
    # cₖ = -- - ----; k = 1, ..., n - 1
    #      4k   4k+2
    c = _chebyshev_c(b, n)

    F_0 = c'U_0[2:n]
    F_s = c'U_s[2:n]

    # ~ Step 6: Compute the value of Dᵅf(s)
    return f(zero(K)) * s^(-α) + ((df_s / (1 - α)) - F_s + F_0) * s^(1 - α)
end

function _chebyshev_method(α::T) where {T<:Number}
    if 0 < α < 1
        return ChebyshevSeriesUnit()
    else
        error(
            """
            No algorithms for α = $(α). Options are:
            - ChebyshevSeriesUnit: 0 < α < 1
            """
        )
    end
end

struct ChebyshevSeries <: ChebyshevDiff end

function fracdiff(f::Union{Function,K}, α::K, s::K, n::Integer, ::ChebyshevSeries) where {K<:Number}
    return fracdiff(f, α, s, n, _chebyshev_method(α))
end