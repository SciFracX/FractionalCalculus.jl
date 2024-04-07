using BenchmarkTools: BenchmarkTools, BenchmarkGroup, @btime, @benchmarkable
using FractionalCalculus
using Statistics: median

const SUITE = BenchmarkGroup()

testf(x) = x^2

SUITE["Caputo"]["CaputoDiethelm"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, CaputoDiethelm())
SUITE["Caputo"]["CaputoTrap"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, CaputoTrap())
SUITE["Caputo"]["CaputoL1"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, CaputoL1())

SUITE["RiemannLiouville"]["RLDiffL1"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, RLDiffL1())
SUITE["RiemannLiouville"]["RLD"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, RLD())

SUITE["GrunwaldLetnikov"]["GLDirect"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, GLDirect())
SUITE["GrunwaldLetnikov"]["GLFiniteDifference"] = @benchmarkable fracdiff(testf, 0.5, 1, 0.01, GLFiniteDifference())

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)

BenchmarkTools.save(joinpath(@__DIR__, "benchmark_results.json"), median(results))