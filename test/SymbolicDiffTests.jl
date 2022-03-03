using FractionalCalculus
using SymbolicUtils
using Test

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
    @eqtest semidiff(x^2) == gamma(3)/gamma(2.5)*(x^1.5)
    @eqtest semidiff(log(x)) == log(4*x) / sqrt(π*x)
    @eqtest semidiff(sqrt(1+x)) == 1/sqrt(pi*x)+atan(sqrt(x))/sqrt(pi)
    @eqtest semidiff(sqrt(1-x)) == 1/sqrt(pi*x)-atanh(sqrt(x))/sqrt(pi)
    @eqtest semidiff(1/sqrt(1+x)) == 1/(sqrt(pi*x)*(1+x))
    @eqtest semidiff(1/sqrt(1-x)) == 1/(sqrt(pi*x)*(1-x))
    @eqtest semidiff(1/(1+x)) == (sqrt(1+x)-sqrt(x)*asinh(sqrt(x)))/(sqrt(pi*x)*(1+x)^(3/2))
    @eqtest semidiff(1/(1-x)) == (sqrt(1-x)+sqrt(x)*asinh(sqrt(x)))/(sqrt(pi*x)*(1-x)^(3/2))
    @eqtest semidiff(1/(1+x)^(3/2)) == (1-x)/(sqrt(pi*x)*(1+x)^2)
    @eqtest semidiff(1/(1-x)^(3/2)) == (1+x)/(sqrt(pi*x)*(1-x)^2)
end