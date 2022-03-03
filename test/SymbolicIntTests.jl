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

@testset "Test INT CONSTANT AND POWERS" begin
    @syms x
    @eqtest semiint(x^3) == 0.51583047638652*(x^3.5)
    @eqtest semiint(log(x)) == 2*log(4*x - 2)*sqrt(0.3183098861837907*x)
    @eqtest semiint(sqrt(1+x)) == sqrt(x/pi) + ((1+x)*arctan(sqrt(x)))/sqrt(pi)
    @eqtest semiint(sqrt(1+x)) == sqrt(x/pi) + ((1-x)*arctanh(sqrt(x)))/sqrt(pi)
    @eqtest semiint(1/sqrt(1+x)) == 2/sqrt(pi)*atan(sqrt(x))
    @eqtest semiint(1/sqrt(1-x)) == 2/sqrt(pi)*atan(sqrt(x))
end