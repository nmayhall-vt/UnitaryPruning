using UnitaryPruning
using Test

@testset "PauliStrings" begin
    
    p1 = PauliString("XIYIY")
    
    p2 = PauliString("XIYZY")

    @test commute(p1,p2)
    
    p1 = PauliString("XIYIX")
    p2 = PauliString("XIYZY")

    @test (commute(p1,p2) == false)
    
    p1 = PauliString("ZIYIX")
    p2 = PauliString("XIYZY")

    @test (commute(p1,p2) == true)
end
