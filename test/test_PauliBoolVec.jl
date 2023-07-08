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

    @test commute(a,b) == false
    
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
    
    @test commute(a,b) == true 

    display(UnitaryPruning.is_diagonal(a))

    
    a = PauliBoolVec("IIZIZIIII")
    v = Vector{Bool}([1,1,1,0,0,0,0,0,0])
    @test expectation_value_sign(a,v) == -1
    
    v = Vector{Bool}([1,1,0,1,0,0,0,0,0])
    @test expectation_value_sign(a,v) == 1
    
    v = Vector{Bool}([1,1,1,1,1,0,0,0,0])
    @test expectation_value_sign(a,v) == 1
end