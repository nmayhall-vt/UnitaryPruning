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




    g = PauliBoolVec(8,X=[1,2], Y=[3], Z=[5,6])
    o = PauliBoolVec(8,X=[3], Y=[1,2], Z=[5,8])

    gmat = to_matrix(g)
    omat = to_matrix(o)
    Imat = to_matrix(PauliBoolVec(8))
    println()
    println("multiplication:")
    e1 = norm(omat*gmat - to_matrix(multiply(o,g)))
    e2 = norm(gmat*gmat - Imat)
    e3 = norm(to_matrix(multiply(g,g)) - Imat)

    @test abs(e1) < 1e-12
    @test abs(e2) < 1e-12
    @test abs(e3) < 1e-12
end