using FractionalCalculus
using Plots

x=collect(0.1:0.1:5)
y=collect(0.1:0.1:5)
z=zeros(50, 50)

for j in range(1, 50, step=1)
    z[j, :] = fracdiff(sin, 0.02*j, y, 0.01, RLDiffApprox())
end

plot(x, y, z, st=:surface)

savefig("./3dexample.png")