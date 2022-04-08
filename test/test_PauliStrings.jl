using UnitaryPruning
using Test

up = UnitaryPruning

@testset "PauliStrings" begin
    
    p1 = PauliString("XIYIY")
    
    p2 = PauliString("XIYZY")

    @test up.commute(p1,p2)
    
    p1 = PauliString("XIYIX")
    p2 = PauliString("XIYZY")
    p3 = PauliString("IIIZZ")
    phase = 1im

    @test (up.commute(p1,p2) == false)
    @test (up.commutator(p1,p2)[1] == phase)
    @test (up.commutator(p1,p2)[2] == p3)
    
    @test (up.commutator(p1,p2)[1] == -up.commutator(p2,p1)[1])
    @test (up.commutator(p1,p2)[2] == up.commutator(p2,p1)[2])
    
    p1 = PauliString("ZIYIX")
    p2 = PauliString("XIYZY")
   
    @test (up.commute(p1,p2) == true)
    
    p1 = PauliString("ZIYIXXYXI")
    p2 = PauliString("XIZIYYIZY")
     
    p3 = PauliString("YIXIZZYYY")
    phase = -1im

    @test (up.commute(p1,p2) == false)
    @test (up.commutator(p1,p2)[1] == phase)
    @test (up.commutator(p1,p2)[2] == p3)


    p = PauliString("IIZIZZIIZZIIIZ")
    @test(up.is_diagonal(p) == true)
    p = PauliString("IIZIZZIIXZIIIZ")
    @test(up.is_diagonal(p) == false)

    p1 = PauliString("YZYIXX")
    p2 = PauliString("XIYZYY")
    p3 = PauliString("ZZIZZZ")
    phase = 1im
    @test (up.commutator(p1,p2)[1] == phase)
    @test (up.commutator(p1,p2)[2] == p3)

    p1 = PauliString("YZYIZX")
    p2 = PauliString("XIYZYY")
    p3 = PauliString("ZZIZXZ")
    phase = -1im
    @test (up.commutator(p1,p2)[1] == phase)
    @test (up.commutator(p1,p2)[2] == p3)


    p3 = PauliString("IIIZZ")
end
