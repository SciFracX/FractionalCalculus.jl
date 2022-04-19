using FractionalCalculus
using Plots, LaTeXStrings

s="\$(x^4)'''=24x\$"

tspan=collect(0:0.01:6)

result=fracdiff(cos, 3, 6, 0.01, RLDiffMatrix())

plot(tspan, result, title=s, legend=:bottomright)