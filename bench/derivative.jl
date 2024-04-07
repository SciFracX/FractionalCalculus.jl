function benchmark_caputo(tag::String, end_tag::String, step_size::Real)
    SUITE["fracdiff"]["cpu"][tag][end_tag] = @benchmarkable fracdiff(x->x^2, 0.5, 1, $step_size, CaputoDiethelm())
    SUITE["fracdiff"]["cpu"][tag][end_tag] = @benchmarkable fracdiff(x->x^2, 0.5, 1, $step_size, CaputoL1())
end

function benchmark_RL(tag::String, end_tag::String, step_size::Real)
    SUITE["fracdiff"]["cpu"][tag][end_tag] = @benchmarkable fracdiff(x->x^2, 0.5, 1, $step_size, RLD())
    SUITE["fracdiff"]["cpu"][tag][end_tag] = @benchmarkable fracdiff(x->x^2, 0.5, 1, $step_size, RLDiffL1())
end