using FractionalCalculus
using Plots, LaTeXStrings

s="\$(x^2)'=2x\$"

tspan=collect(0:0.01:3)

result=fracdiff(x->x^2, 1, 3, 0.01, RLDiff_Matrix())

plot(tspan, result, title=s, legend=:bottomright)