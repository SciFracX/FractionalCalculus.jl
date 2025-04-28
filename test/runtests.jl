using FractionalCalculus
using SymbolicUtils
using SpecialFunctions
using Test

@testset "FractionalCalculus" begin
    include("aqua.jl")
    include("./Derivative.jl")
    include("./Integral.jl")
    include("./SymbolicDiffTests.jl")
    include("./SymbolicIntTests.jl")
end