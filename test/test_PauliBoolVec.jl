using UnitaryPruning
using Test
using BenchmarkTools 

@testset "PauliMask" begin

    a = PauliBoolVec("XYZIXY")
    b = PauliBoolVec("YXXYZZ")
    c = a*b
    cref = PauliBoolVec("ZZYYYX")
    cref.θ = 2
    display(c)
    display(cref)
    @test c == cref 

end