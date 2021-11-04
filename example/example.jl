using LinearAlgebra, InvertedIndices
using Plots


"""
Generating elements in Matrix.
"""
function omega(n, p)
    omega = zeros(n+1)

    omega[1]=1
    for i in range(1, n, step=1)
        omega[i+1]=(1-(p+1)/i)*omega[i]
    end
    
    return omega

end

function B(N, p, h)
    result=zeros(N, N)
    temp=omega(N, p)

    for i in range(1, N, step=1)
        result[i, 1:i]=reverse(temp[1:i])
    end

    return h^(-p)*result
end

tspan=collect(0:0.001:1)
result=B(1001, 0.5, 0.001)*tspan.^2
#print(result[end])
#Plots.plot(tspan, result)