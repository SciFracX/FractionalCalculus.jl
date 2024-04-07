using BenchmarkTools: BenchmarkTools, BenchmarkGroup, @btime, @benchmarkable
using FractionalCalculus
using Statistics: median

const SUITE = BenchmarkGroup()

testf(x) = x^2
for h in (0.01, 0.001)
    SUITE["Caputo"]["CaputoDiethelm"] = @benchmarkable fracdiff(testf, 0.5, 1, h, CaputoDiethelm())
    SUITE["Caputo"]["CaputoTrap"] = @benchmarkable fracdiff(testf, 0.5, 1, h, CaputoTrap())
    SUITE["Caputo"]["CaputoL1"] = @benchmarkable fracdiff(testf, 0.5, 1, h, CaputoL1())

    SUITE["RiemannLiouville"]["RLDiffL1"] = @benchmarkable fracdiff(testf, 0.5, 1, h, RLDiffL1())
    SUITE["RiemannLiouville"]["RLD"] = @benchmarkable fracdiff(testf, 0.5, 1, h, RLD())

    SUITE["GrunwaldLetnikov"]["GLFiniteDifference"] = @benchmarkable fracdiff(testf, 0.5, 1, h, GLFiniteDifference())
end

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)

BenchmarkTools.save(joinpath(@__DIR__, "benchmark_results.json"), median(results))