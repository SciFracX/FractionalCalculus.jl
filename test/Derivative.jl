using FractionalCalculus
using Test

@testset "Test Caputo Fractional Derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, 0, 1, 1e-8, Caputo())[1], 2/sqrt(pi); atol = 1e-5)
    @test isapprox(fracdiff(x->x^5, 0.5, 0, 3.2, 1e-8, Caputo())[1], 4.300306216488329e2; atol = 1e-5)
    @test isapprox(fracdiff(x->x^5, 0.25, 0, 3.2, 1e-8, Caputo())[1], 382.1228113951680; atol = 1e-5)
    @test isapprox(fracdiff(x->(x-3)^3+(x+2)^2+5, 0.5, 0, 0.83, 1e-8, Caputo())[1], 23.899930581118699; atol = 1e-5)
    @test isapprox(fracdiff(x->(x-3)^3+(x+2)^2+5, 0.25, 0, 0.83, 1e-8, Caputo())[1], 22.96354626572943; atol = 1e-5)
    @test isapprox(fracdiff(x->cos(x), 0.5, 0, 0.34, 1e-8, Caputo())[1], -0.147174774607683; atol = 1e-5)
    @test isapprox(fracdiff(x->cos(x), 0.25, 0, 0.34, 1e-8, Caputo())[1], -0.093074344902431; atol = 1e-5)
    @test isapprox(fracdiff(x->exp(-x^2), 0.5, 0, 2.3, 1e-8, Caputo())[1], -0.505891017154289; atol = 1e-5)
    @test isapprox(fracdiff(x->exp(-x^2), 0.25, 0, 2.3, 1e-8, Caputo())[1], -0.762927553656252; atol = 1e-5)
end

@testset "Test Grunwald-Letnikov Fractional Derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, 0, 0.5, GL())[1], 0.7978845583186518; atol = 1e-5)
end