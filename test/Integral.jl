using FractionalCalculus
using Test

@testset "Test Direct Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 0, 1.23, 1e-8, RL_Direct())[1], 1.163931646862977; atol = 0.00001)
    @test isapprox(fracint(x->x^4, 0.25, 0, 1.23, 1e-8, RL_Direct())[1], 1.64294134724409; atol = 0.00001)
    @test isapprox(fracint(x->x^4, 0.15, 0, 1.23, 1e-8, RL_Direct())[1], 1.878978485456182; atol = 0.00001)
    @test isapprox(fracint(x->sin(x), 0.5, 0, 0.97, 1e-8, RL_Direct())[1], 0.6442601778843754; atol = 0.00001)
    @test isapprox(fracint(x->sin(x), 0.25, 0, 0.97, 1e-8, RL_Direct())[1], 0.7448298999485293; atol = 0.00001)
    @test isapprox(fracint(x->sin(x), 0.15, 0, 0.97, 1e-8, RL_Direct())[1], 0.7802555196986999; atol = 0.00001)
end

@testset "Test Piecewise Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-4, RL_Piecewise()), 1.163931646862977; atol = 1e-5)
    @test isapprox(fracint(x->x^4, 0.25, 1.23, 1e-4, RL_Piecewise()), 1.64294134724409; atol = 1e-5)
    @test isapprox(fracint(sin, 0.5, 0.97, 1e-4, RL_Piecewise()), 0.6442601778843754; atol = 1e-5)
end

@testset "Test Simpson Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-7, RLInt_Simpson()), 1.163931646862977; atol = 1e-3)
end

@testset "Test Trapezoidal Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-4, RLInt_Trapezoidal()), 1.163931646862977; atol = 1e-4)
end

@testset "Test Rectangular Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-7, RLInt_Rectangular()), 1.163931646862977; atol = 1e-4)
end