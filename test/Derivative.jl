using FractionalCalculus, SpecialFunctions
using Test

half_derivative_of_sqrt_x_at_1 = sqrt(π)/2
half_derivative_of_x_square_minus_1_at_1 = 5/(3*sqrt(π))
half_derivative_of_exp_at_1 = ℯ*erf(1)+1/sqrt(π)

@testset "Test Caputo Direct Fractional Derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, 0, 1, 1e-4, CaputoDirect())[1], 2/sqrt(pi); atol = 1e-5)
    @test isapprox(fracdiff(x->x^5, 0.5, 0, 3.2, 1e-4, CaputoDirect())[1], 4.300306216488329e2; atol = 1e-5)
    @test isapprox(fracdiff(x->x^5, 0.25, 0, 3.2, 1e-4, CaputoDirect())[1], 382.1228113951680; atol = 1e-5)
    @test isapprox(fracdiff(x->(x-3)^3+(x+2)^2+5, 0.5, 0, 0.83, 1e-4, CaputoDirect())[1], 23.899930581118699; atol = 1e-5)
    @test isapprox(fracdiff(x->(x-3)^3+(x+2)^2+5, 0.25, 0, 0.83, 1e-4, CaputoDirect())[1], 22.96354626572943; atol = 1e-5)
    @test isapprox(fracdiff(x->cos(x), 0.5, 0, 0.34, 1e-4, CaputoDirect())[1], -0.147174774607683; atol = 1e-5)
    @test isapprox(fracdiff(x->cos(x), 0.25, 0, 0.34, 1e-4, CaputoDirect())[1], -0.093074344902431; atol = 1e-5)
    @test isapprox(fracdiff(x->exp(-x^2), 0.5, 0, 2.3, 1e-4, CaputoDirect())[1], -0.505891017154289; atol = 1e-5)
    @test isapprox(fracdiff(x->exp(-x^2), 0.25, 0, 2.3, 1e-4, CaputoDirect())[1], -0.762927553656252; atol = 1e-5)

    @test isapprox(fracdiff(x->x, 0.5, 0, [1, 2, 3], 0.0001, CaputoDirect()), [2/sqrt(pi), 1.5957691216057408, 1.9544100476116828]; atol=1e-3)
end

@testset "Test Caputo Trapezoidal Fractional Derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 1e-4, CaputoTrap()), 2/sqrt(pi); atol = 1e-5)
    @test isapprox(fracdiff(x->x^5, 0.5, 3.2, 1e-4, CaputoTrap()), 4.300306216488329e2; atol = 1e-5)
    @test isapprox(fracdiff(cos, 0.5, 0.34, 1e-4, CaputoTrap()), -0.147174774607683; atol = 1e-5)
    @test isapprox(fracdiff(x->exp(-x^2), 0.25, 2.3, 1e-3, CaputoTrap()), -0.762927553656252; atol = 1e-5)
end

@testset "Test Caputo Piecewise Fractional Derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, [1, 2, 3], 1e-4, CaputoTrap()), [1.1283791670557781, 1.5957691213491374, 1.9544100474864603])
end

@testset "Test Caputo Diethelm Fractional Derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 1e-4, CaputoDiethelm()), 2/sqrt(pi); atol = 1e-5)
    @test isapprox(fracdiff(x->x^5, 0.5, 3.2, 1e-5, CaputoDiethelm()), 4.300306216488329e2; atol = 1e-5)
    @test isapprox(fracdiff(x->x, 0.5, [1, 2, 3], 0.0001, CaputoDiethelm()), [1.1283791670955023, 1.5957691216057386, 1.954410047611745])
end

@testset "Test Caputo High Order Fractional Derivative" begin
    @test isapprox(fracdiff(x->x^8, 0.5, 1, 0.05, 3, CaputoHighOrder()), 2.8542562386883445; atol=1e-4)
end

@testset "Test Grunwald-Letnikov Fractional Derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, 0, 0.5, GLDirect())[1], 0.7978845583186518; atol = 1e-5)
    @test isapprox(fracdiff(x->x^5, 0.5, 0, 3.2, GLDirect())[1], 4.300306216488329e2; atol = 1e-5)
    @test isapprox(fracdiff(x->x, 0.5, collect(0:0.01:1), GLHighPrecision())[end], 2/sqrt(pi); atol = 1e-4)
end

@testset "Test GL_Multiplicative_Additive()" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.009, GLMultiplicativeAdditive()), 2/sqrt(pi); atol=1e-2)
end

@testset "Test GL_Lagrange_Three_Point_Interp" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.009, GLLagrangeThreePointInterp()), 2/sqrt(pi); atol=1e-2)
end

@testset "Test GL_Finite_Difference" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.009, GLFiniteDifference()), 2/sqrt(pi); atol=1e-2)
end

@testset "Test Grunwald Letnikov Lagrange_Three_Point_Interp" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.01, GLLagrangeThreePointInterp()), 2/sqrt(pi); atol=1e-4)
end

@testset "Test RLDiff Matrix Fractional Deriavtive" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.001, RLDiffMatrix())[end], 2/sqrt(pi); atol = 1e-3)
end

@testset "Test RLDiffL1" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.001, RLDiffL1()), 2/sqrt(pi); atol = 1e-4)

    # Test Array input end point
    result = fracdiff(x->x, 0.5, [1, 2, 3], 0.001, RLDiffL1())
    @test isapprox(result[1], 2/sqrt(pi); atol=1e-4)
    @test isapprox(result[2], 1.5957691216057306; atol=1e-4)
    @test isapprox(result[3], 1.954410047611678; atol=1e-4)
end

@testset "Test RL Linear Spline Interpolation" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.0001, RLLinearSplineInterp()), 2/sqrt(pi); atol = 1e-4)

    # Test Array input end point
    result = fracdiff(x->x, 0.5, [1, 2, 3], 0.001, RLLinearSplineInterp())
    @test isapprox(result[1], 2/sqrt(pi); atol=1e-4)
    @test isapprox(result[2], 1.5957691216057306; atol=1e-4)
    @test isapprox(result[3], 1.954410047611678; atol=1e-4)
end


@testset "Test RL_G1 method" begin
    @test isapprox(fracdiff(x->x, 0.5, 0, 1, 0.006, RLG1()), 2/sqrt(pi); atol=1e-2)
end

@testset "Test RL_D method" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.001, RLD()), 2/sqrt(pi); atol=1e-4)
end

@testset "Test Hadamard Fractional Derivative" begin
    @test isapprox(fracdiff(log, 0.5, 1, 2, 1/4000000, HadamardLRect()), 0.9391460; atol=1e-4)
    @test isapprox(fracdiff(log, 0.3, 1, 2, 1/4000000, HadamardLRect()), 0.8514935; atol=1e-4)
    @test isapprox(fracdiff(log, 0.5, 1, 2, 1/4000000, HadamardRRect()), 0.9391460; atol=1e-4)
    @test isapprox(fracdiff(log, 0.3, 1, 2, 1/4000000, HadamardRRect()), 0.8514936; atol=1e-4)
    @test isapprox(fracdiff(log, 0.5, 1, 2, 1/4000000, HadamardTrap()), 0.9391460; atol=1e-4)
    @test isapprox(fracdiff(log, 0.3, 1, 2, 1/4000000, HadamardTrap()), 0.8514935; atol=1e-4)
end

@testset "Test Riesz Symmetric fractional derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.5, RieszSymmetric()), [0.2651650429449553; 0.2651650429449553; -0.39774756441743303]; atol=1e-2)
end

@testset "Test Caputo-Fabrizio fractional derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.00001, CaputoFabrizioAS()), 0.9887564512257243; atol=1e-4)
    @test isapprox(fracdiff(x->x^2, 0.5, 1, 0.0000001, CaputoFabrizioAS()), 1.1508665775147289; atol=1e-6)
end

@testset "Test Atangana-Baleanu Caputo sense fractional derivative" begin
    @test isapprox(fracdiff(x->x, 0.5, 1, 0.00001, AtanganaSeda()), -0.8696378200415389; atol=1e-3)
end

@testset "Test Chebyshev sense fractional derivative" begin
    let f(x) = x^7, α = 0.4, n = 10
        @test_broken isapprox(fracdiff(f, 0.3, α, n, ChebyshevSeries()), 0.00509416; atol=1e-5)
        @test_broken isapprox(fracdiff(f, 0.5, α, n, ChebyshevSeries()), 0.01236690; atol=1e-5)
        @test_broken isapprox(fracdiff(f, 0.7, α, n, ChebyshevSeries()), 0.03689920; atol=1e-5)
    end

    let f(x) = x^7, α = 0.4, β = 0.7, n = 10
        @test_broken isapprox(fracdiff(f, 0.3, α, n, ChebyshevSeries()), 0.02895030; atol=1e-5)
        @test_broken isapprox(fracdiff(f, 0.5, α, n, ChebyshevSeries()), 0.06579140; atol=1e-5)
        @test_broken isapprox(fracdiff(f, 0.5, β, n, ChebyshevSeries()), 0.81628300; atol=1e-5)
        @test_broken isapprox(fracdiff(f, 0.7, α, n, ChebyshevSeries()), 0.18334300; atol=1e-5)
    end
end

@testset "Test macros" begin
    @test isapprox(@fracdiff(x->x, 0.5, 1), 1.1283791670955126; atol=1e-4)
    @test isapprox(@fracdiff(x->x^5, 0.5, 3.2), 4.300306216488329e2; atol=1e-2)
end


@testset "Test auxiliary functions" begin
    @test isapprox(RieszMatrix(0.5, 2, 0.1), [-1.58114 1.3835; 1.3835 -1.58114]; atol=1e-3)
    @test isapprox(B(2, 2), [1.0 0.0; -2.0 1.0]; atol=1e-3)
    @test isapprox(genfun(2), [1.5 -2 0.5]; atol=1e-3)
end

@testset "Test Mittag Leffler function" begin
    @test isapprox(mittleff(2, 2, 1), sinh(1); atol=1e-5)
end