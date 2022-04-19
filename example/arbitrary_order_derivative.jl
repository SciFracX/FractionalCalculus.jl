using FractionalCalculus
using Plots, LaTeXStrings

target = gamma(6)/gamma(2.4)*tspan.^1.4


s="\$D^{3.6}x^5\$"

tspan=collect(0:0.01:6)

result=fracdiff(x->x^5, 3.6, 6, 0.01, RLDiffMatrix())

plot(tspan, result, title=s, legend=:bottomright, label="Numerical")
plot!(tspan, target, lw=3, ls=:dash, label="Analytical")
savefig("./arbitrary_order_derivative.png")