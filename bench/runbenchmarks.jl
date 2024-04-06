using BenchmarkTools: BenchmarkTools, BenchmarkGroup, @btime, @benchmarkable
using FractionalCalculus
using Statistics: median

const SUITE = BenchmarkGroup()

testf(x) = x^2

SUITE["Caputo"]["CaputoDiethelm"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, CaputoDiethelm())
SUITE["Caputo"]["CaputoTrap"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, CaputoTrap())
SUITE["Caputo"]["CaputoL1"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, CaputoL1())

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)

BenchmarkTools.save(joinpath(@__DIR__, "benchmark_results.json"), median(results))