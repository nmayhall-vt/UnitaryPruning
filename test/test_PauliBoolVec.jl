using UnitaryPruning
using Test
using BenchmarkTools 

@testset "PauliBoolVec" begin

    a = PauliBoolVec("XYZIXY")
    b = PauliBoolVec("YXXYZZ")
    c = PauliBoolVec("ZZYYYX")
    c.θ = 2
    display(c)
    display(a*b)
    @test c == a*b 

    @test commute_check(a,b) == false
    
    c.θ = 0
    @test c == b*a 

    a = PauliBoolVec("XYYZIZYZIIXYIIYZI")
    b = PauliBoolVec("YYYIIZIIZIZYXIXYZ")
    c = PauliBoolVec("ZIIZIIYZZIYIXIZXZ")
    c.θ = 0
    display(a)
    display(b)
    display(c)
    display(a*b)
    @test c == a*b 
    @test c == b*a
    
    @test commute_check(a,b) == true 

end