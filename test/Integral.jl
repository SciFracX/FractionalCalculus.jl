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

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-4, RL_Piecewise())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
end

@testset "Test RLInt_Approx" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-4, RLInt_Approx()), 1.163931646862977; atol = 1e-4)

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-4, RLInt_Approx())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
end

@testset "Test RLInt_Matrix" begin
    @test isapprox(fracint(x->x, 0.5, 1, 0.5, RLInt_Matrix()), [0; 0.3535533905932738; 0.8838834764831844]; atol=1e-3)
end

@testset "Test RL_LinearInterp" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-4, RL_LinearInterp()), 1.163931646862977; atol = 1e-4)

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-4, RL_LinearInterp())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
end

@testset "Test Simpson Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-7, RLInt_Simpson()), 1.163931646862977; atol = 1e-2)

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-7, RLInt_Simpson())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-2)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-2)
    @test isapprox(result[3], 64.3280200256207; atol = 1)
end

@testset "Test Trapezoidal Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-7, RLInt_Trapezoidal()), 1.163931646862977; atol = 1e-4)
#=
    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-7, RLInt_Trapezoidal())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
    =#
end

@testset "Test Rectangular Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-6, RLInt_Rectangular()), 1.163931646862977; atol = 1e-3)

    @test isapprox(fracint(x->x^4, 0.5, [1.23, 2, 3], 0.001, RLInt_Rectangular()), [1.1617493910169032; 10.363130001411157; 64.27898282619807]; atol=1e-2)
end

@testset "Test Cubic Spline Interpolation Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 0.001, RLInt_Cubic_Spline_Interp()), 1.163931646862977; atol = 1e-5)
    @test isapprox(fracint(sin, 0.5, 0.97, 0.01, RLInt_Cubic_Spline_Interp()), 0.6442601778843754; atol = 1e-6)

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 0.01, RLInt_Cubic_Spline_Interp())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
end

@testset "Test fractional integral macros" begin
    @test isapprox(@fracint(x->x^4, 0.5, 1.23), 1.163931646862977; atol = 1e-4)
    @test isapprox(@semifracint(x->x^4, 1.23), 1.163931646862977; atol = 1e-4)
end