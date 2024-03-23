using FractionalCalculus
using Test

@testset "Test Direct Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 0, 1.23, 1e-4, RLDirect())[1], 1.163931646862977; atol = 0.00001)
    @test isapprox(fracint(x->x^4, 0.25, 0, 1.23, 1e-4, RLDirect())[1], 1.64294134724409; atol = 0.00001)
    @test isapprox(fracint(x->x^4, 0.15, 0, 1.23, 1e-4, RLDirect())[1], 1.878978485456182; atol = 0.00001)
    @test isapprox(fracint(x->sin(x), 0.5, 0, 0.97, 1e-4, RLDirect())[1], 0.6442601778843754; atol = 0.00001)
    @test isapprox(fracint(x->sin(x), 0.25, 0, 0.97, 1e-4, RLDirect())[1], 0.7448298999485293; atol = 0.00001)
    @test isapprox(fracint(x->sin(x), 0.15, 0, 0.97, 1e-4, RLDirect())[1], 0.7802555196986999; atol = 0.00001)
end

@testset "Test Piecewise Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-4, RLPiecewise()), 1.163931646862977; atol = 1e-5)
    @test isapprox(fracint(x->x^4, 0.25, 1.23, 1e-4, RLPiecewise()), 1.64294134724409; atol = 1e-5)
    @test isapprox(fracint(sin, 0.5, 0.97, 1e-4, RLPiecewise()), 0.6442601778843754; atol = 1e-5)

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-4, RLPiecewise())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
end

@testset "Test RLInt_Approx" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-4, RLIntApprox()), 1.163931646862977; atol = 1e-4)

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-4, RLIntApprox())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
end

@testset "Test RLInt_Matrix" begin
    @test isapprox(fracint(x->x, 0.5, 1, 0.5, RLIntMatrix()), [0; 0.3535533905932738; 0.8838834764831844]; atol=1e-3)
end

@testset "Test RL_LinearInterp" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-4, RLLinearInterp()), 1.163931646862977; atol = 1e-4)

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-4, RLLinearInterp())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
end

@testset "Test Simpson Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-7, RLIntSimpson()), 1.163931646862977; atol = 1e-2)

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-7, RLIntSimpson())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-2)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-2)
    @test isapprox(result[3], 64.3280200256207; atol = 1)
end

@testset "Test Trapezoidal Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-7, RLIntTrapezoidal()), 1.163931646862977; atol = 1e-4)
#=
    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 1e-7, RLInt_Trapezoidal())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
    =#
end

@testset "Test Rectangular Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 1e-6, RLIntRectangular()), 1.163931646862977; atol = 1e-3)

    @test isapprox(fracint(x->x^4, 0.5, [1.23, 2, 3], 0.001, RLIntRectangular()), [1.1617493910169032; 10.363130001411157; 64.27898282619807]; atol=1e-2)
end

@testset "Test Cubic Spline Interpolation Fractional Integral" begin
    @test isapprox(fracint(x->x^4, 0.5, 1.23, 0.001, RLIntCubicSplineInterp()), 1.163931646862977; atol = 1e-5)
    @test isapprox(fracint(sin, 0.5, 0.97, 0.01, RLIntCubicSplineInterp()), 0.6442601778843754; atol = 1e-6)

    result = fracint(x->x^4, 0.5, [1.23, 2, 3], 0.01, RLIntCubicSplineInterp())
    @test isapprox(result[1], 1.163931646862977; atol = 1e-3)
    @test isapprox(result[2], 10.375035632494745; atol = 1e-3)
    @test isapprox(result[3], 64.3280200256207; atol = 1e-3)
end

@testset "Test fractional integral macros" begin
    @test isapprox(@fracint(x->x^4, 0.5, 1.23), 1.163931646862977; atol = 1e-4)
    @test isapprox(@semifracint(x->x^4, 1.23), 1.163931646862977; atol = 1e-4)
end



@testset "Test Hadamard integral" begin
    (a, b)=fracint(sin, 0.5, 1, 2, 1/3, 0.1, HadamardMat())

    @test isapprox(a, [ 1.0
    1.023136962522694
    1.077552576065
    1.1626221795557758
    1.2764357600635179
    1.414213562373095
    1.5668630279525204
    1.7202493081322234
    1.85605792647593
    1.954772501883529
    2.0]; atol=1e-3)
    @test isapprox(b, [Inf
    3.2545536740283265
    1.9422879792123904
    1.5034729488252023
    1.291625449283718
    1.1492602172225017
    1.0157975575699352
    0.8713525896641822
    0.725484956837982
    0.6074878993825767
    0.5500616274393197]; atol=1e-3)
end