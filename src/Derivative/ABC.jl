abstract type AtanganaBaleanu <: FracDiffAlg end

struct AtanganaSeda <: AtanganaBaleanu end

function fracdiff(f, α, point, h, ::AtanganaSeda)
    n::Int = round(Int, point/h)
    result = zero(Float64)
    AB=1-α+α/gamma(α)
    for j in 0:n
        result += (f((j+1)*h)-f(j*h))/h*((n-j)*h*mittleff(α, 2, -α/(1-α)*(n-j)^α*(h)^α)-(n-j+1)*h*mittleff(α, 2, -α/(1-α)*(n-j+1)^α*h^α))
    end
    return result*AB/(1-α)
end