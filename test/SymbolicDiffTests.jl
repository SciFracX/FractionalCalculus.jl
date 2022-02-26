using FractionalCalculus
using SymbolicUtils
using Test

# == / != syntax is nice, let's use it in tests
macro eqtest(expr)
    @assert expr.head == :call && expr.args[1] in [:(==), :(!=)]
    if expr.args[1] == :(==)
        :(@test isequal($(expr.args[2]), $(expr.args[3])))
    else
        :(@test !isequal($(expr.args[2]), $(expr.args[3])))
    end |> esc
end

@testset "Test DIFF CONSTANT AND POWERS" begin
    @syms x
    @eqtest semidiff(x^2) == 1.5045055561273502*(x^1.5)
    @eqtest semidiff(log(x)) == log(4x) / sqrt(π*x)
end