using Aqua
using Test
using FractionalCalculus

@testset "Aqua" begin
    Aqua.test_all(FractionalCalculus)
end