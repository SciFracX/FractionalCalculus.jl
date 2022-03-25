abstract type HadamardInt <: FracIntAlg end

"""
    fracint(f, α, a, b, mu, N, HadamardMat())

Compute Hadamard fractional integral.

### References

```tex
@article{Xu2019SpectralCM,
  title={Spectral Collocation Method for Fractional Differential/Integral Equations with Generalized Fractional Operator},
  author={Qinwu Xu and Zhoushun Zheng},
  journal={International Journal of Differential Equations},
  year={2019}
}
```
"""
struct HadamardMat <: HadamardInt end

function fracint(g, α, a, b, mu, N, ::HadamardMat)
    w(x)=x.^mu
    invz(x)=exp.(x)
    (rx, r, D) = GFracMat(log(a), log(b), N, α, w, invz)
    Ga=a^mu.*g.(a)/gamma(1-α)*rx.^(-mu).*(log.(rx./a)).^(-α)
    G=g.(rx)
    H=Ga+D*G
    return rx, H
end

function GFracMat(za, zb, N, α, w, invz)
    r0=((JacobiGL(0, 0, N).+1)./2)
    col_pts=za.+r0*(zb-za)
    col_pt = deepcopy(col_pts)
    if α>0
    FM=Jac_Frac_Diff_Shift(col_pts,0,0,N,α,za,zb,1)
    else
        FM=(2/(zb-za))^α*Jac_Frac_Int(r0*2-1, 0, 0, N, -α, 1)
        FM[1, :] .= 0
    end

    r1 = r0
    interp_pts = r1.*2 .-1
    Vand = Jac_Va(interp_pts, 0, 0, N, 1)

    rx=invz.(col_pt)
    WL = diagm(1 ./w.(rx))
    WL[1, 1]=0
    WR = diagm(w.(rx))
    M=WL*(FM/Vand)*WR
    return rx, col_pt, M
end

function JacobiGL(α,β,N)
    x = zeros(N+1, 1)
    if N==1
        x[1] = -1.0
        x[2] = 1.0
        return x
    end
    
    (xint, _) = JacobiGQ(α+1,β+1,N-2)
    x = [-1 xint' 1]'
    return x
end
    
function JacobiGQ(α, β, N)
    if N==0
        x[1] = -(α-β)/(α+β+2)
        w[1] = 2
        return x, w
    end

    J = zeros(N+1)
    h1 = 2*collect(0:N).+α.+β
    J = diagm(-1/2*(α^2-β^2)./(h1.+2)./h1) + 
        diagm(1 => 2 ./(h1[1:N].+2).*sqrt.(collect(1:N).*(collect(1:N).+α.+β).*(collect(1:N).+α).*(collect(1:N).+β)./(h1[1:N].+1)./(h1[1:N].+3)))
    if α+β < 10*eps()
        J[1, 1] = 0.0
    end
    J = J + J'
    
    (D, V) = eigen(J)
    x = diagm(D)

    w = V[1, :].^2*2^(α+β+1)/(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(α+β+1)
    return diag(x), w
end

function Jac_Frac_Diff(x, a, b,N, α, output_type)
    if size(x, 2) !== 1
        x=x'
    end
    L=length(x)
    
    k = ceil(Int64, α)
    d = gamma.(collect(0:N).+(k+a+b+1))./(2^k*gamma.(collect(0:N).+(a+b+1)))
    
    if output_type==0
        if N < k
            JFD = zeros(size(x))
            return JFD
        end
        JFD = d[N+1]*Jac_Frac_Int(x,a+k,b+k,N-k,k-α,0)
    else
        JFD = zeros(L, N+1)
        JFD[:, k+1:N+1] = Jac_Frac_Int(x, a+k, b+k, N-k, k-α, 1)
        for j = k+1:N+1
            JFD[:, j] = d[j]*JFD[:, j]
        end
    end
    return JFD
end

function Jac_Frac_Diff_Shift(x,a,b,N,α,min,Max,output_type)
    L=Max-min
    x = 2 .*(x.-min)./L.-1
    FJ = (2/L)^α*Jac_Frac_Diff(x, a, b, N, α, output_type)
    return FJ
end


function Jac_Frac_Int(x, a, b, N::Int64, α, output_type)
    if size(x,2) !== 1
        x=x'
    end
    L=length(x)
    J = zeros(L, Int64(N+1))
    @. J[:, 1] = (x+1)^α/gamma(α+1)
    if N==0
        return J
    end
    @. J[:, 2] = (a+b+2)/2*(x.*(x+1).^α/gamma(α+1)-α*(x+1).^(α+1)/gamma(α+2))+(a-b)/2*J[:, 1]
    if N==1
        return J
    end
    Jac_1 = Jac_Va(-1, a, b, N, 1)
    Aj(a, b, j) = @. (2*j+a+b+1)*(2*j+a+b+2)/(2*(j+1)*(j+a+b+1))
    Bj(a, b, j) = @. (b^2-a^2)*(2*j+a+b+1)/(2*(j+1)*(j+a+b+1)*(2*j+a+b))
    Cj(a, b, j) = @. (j+a)*(j+b)*(2*j+a+b+2)/((j+1)*(j+a+b+1)*(2*j+a+b))
    Ahatj(a, b, j) = @. -2*(j+a)*(j+b)/((j+a+b)*(2*j+a+b)*(2*j+a+b+1))
    Bhatj(a, b, j) = @. 2*(a-b)/((2*j+a+b)*(2*j+a+b+2))
    Chatj(a, b, j) = @. 2*(j+a+b+1)/((2*j+a+b+1)*(2*j+a+b+2))

    for n=2:N
        @. J[:, n+1] = (Aj(a,b,n-1)*x-Bj(a,b,n-1)-α*Aj(a,b,n-1)*Bhatj(a,b,n-1))*J[:, n]-(Cj(a, b, n-1) + α*Aj(a, b, n-1)*Ahatj(a, b, n-1)).*J[:, n-1] + Aj(a, b, n-1)*(Ahatj(a, b, n-1)*Jac_1[n-1]+Bhatj(a,b,n-1)*Jac_1[n]+Chatj(a,b,n-1)*Jac_1[n+1])*(x+1).^α/gamma(α)
        @. J[:, n+1] =J[:, n+1]/(1+α*Aj(a,b,n-1)*Chatj(a,b,n-1))
    end
    if output_type==0
        J=J[:, N+1]
    end
    return J
end

function Jac_Va(x, a, b, N, output_type)
    if size(x,2) !== 1
        x=x'
    end
    L=length(x)
    J = ones(L, N+1)
    if N==0
        J=J[:, 1]
        return J
    end
    @. J[:, 2]=(1+a)+(a+b+2)*(x/2 -1/2)
    if N==1
        if output_type==0
            J=J[:, 2]
            return J
        else
            J=J[:, 1:2]
            return J
        end
    end
    for n=2:N
        @. J[:, n+1] = ((2*n+a+b-1)*((2*n+a+b)*(2*n+a+b-2)*x+a^2-b^2).*J[:, n]-2*(n+a-1)*(n+b-1)*(2*n+a+b)*J[:, n-1])/(2*n*(n+a+b)*(2*n+a+b-2))
    end
    if output_type==0
        J=J[:, N+1]
    end
    return J
end
