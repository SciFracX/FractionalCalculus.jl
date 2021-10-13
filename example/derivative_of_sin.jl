using FractionalCalculus
using Plots

tspan=collect(0:0.01:15)

examplederivative = []

for i in range(0.1, 0.9, step=0.1)
    push!(examplederivative, fracdiff(sin, i, tspan, 0.001, RLDiff_Approx()))
end

plot(tspan, examplederivative, title="Different order of derivative", linewidth=3, label=["α=0.1" "α=0.2" "α=0.3" "α=0.4" "α=0.5" "α=0.6" "α=0.7" "α=0.8" "α=0.9" "α=1"], legend=:bottomright)