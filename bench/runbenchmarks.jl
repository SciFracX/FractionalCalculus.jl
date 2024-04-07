using BenchmarkTools: BenchmarkTools, BenchmarkGroup, @btime, @benchmarkable
using FractionalCalculus
using Statistics: median

const SUITE = BenchmarkGroup()

include("derivative.jl")

for step_size in (0.1, 0.01, 0.001)
    benchmark_caputo("Caputo", "CaputoDiethelm", step_size)
    benchmark_caputo("Caputo", "CaputoL1", step_size)
    benchmark_RL("RiemannLiouville", "RLD", step_size)
    benchmark_RL("RiemannLiouville", "RLDiffL1", step_size)
end

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)

BenchmarkTools.save(joinpath(@__DIR__, "benchmark_results.json"), median(results))