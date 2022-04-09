using UnitaryPruning
using Test
using BenchmarkTools 

@testset "PauliStrings" begin
    p1 = PauliString("XIYIY")
    
    p2 = PauliString("XIYZY")

    @test UnitaryPruning.commute(p1,p2)
    
    p1 = PauliString("XIYIX")
    p2 = PauliString("XIYZY")
    p3 = PauliString("IIIZZ")
    phase = 1im

    @test (UnitaryPruning.commute(p1,p2) == false)
    @test (UnitaryPruning.commutator(p1,p2)[1] == phase)
    @test (UnitaryPruning.commutator(p1,p2)[2] == p3)
    
    @test (UnitaryPruning.commutator(p1,p2)[1] == -UnitaryPruning.commutator(p2,p1)[1])
    @test (UnitaryPruning.commutator(p1,p2)[2] == UnitaryPruning.commutator(p2,p1)[2])
    
    p1 = PauliString("ZIYIX")
    p2 = PauliString("XIYZY")
   
    @test (UnitaryPruning.commute(p1,p2) == true)
    
    p1 = PauliString("ZIYIXXYXI")
    p2 = PauliString("XIZIYYIZY")
     
    p3 = PauliString("YIXIZZYYY")
    phase = -1im

    @test (UnitaryPruning.commute(p1,p2) == false)
    @test (UnitaryPruning.commutator(p1,p2)[1] == phase)
    @test (UnitaryPruning.commutator(p1,p2)[2] == p3)


    p = PauliString("IIZIZZIIZZIIIZ")
    @test(UnitaryPruning.is_diagonal(p) == true)
    p = PauliString("IIZIZZIIXZIIIZ")
    @test(UnitaryPruning.is_diagonal(p) == false)

    p1 = PauliString("YZYIXX")
    p2 = PauliString("XIYZYY")
    p3 = PauliString("ZZIZZZ")
    phase = 1im
    @test (UnitaryPruning.commutator(p1,p2)[1] == phase)
    @test (UnitaryPruning.commutator(p1,p2)[2] == p3)

    p1 = PauliString("YZYIZX")
    p2 = PauliString("XIYZYY")
    p3 = PauliString("ZZIZXZ")
    phase = -1im
    @test (UnitaryPruning.commutator(p1,p2)[1] == phase)
    @test (UnitaryPruning.commutator(p1,p2)[2] == p3)

    # pauli bitstrings
    pbs1 = PauliBitString{UInt16,8}(0,0)
    pbs1 = PauliBitString(p1, T=UInt16)
    pbs2 = PauliBitString(p2, T=UInt16)
    pbs3 = PauliBitString(p3, T=UInt16)

    @test UnitaryPruning.n_non_eye(pbs3) == 5
   
    println(" commute: ", UnitaryPruning.commute(pbs1,pbs2))
    #@btime UnitaryPruning.commute($pbs1, $pbs2)

    # Check that we can interconvert
    pbs1 = PauliBitString{UInt,20}(367,23)
    ps1 = PauliString("XIYZYYIIZXYZ")
    @test pbs1 == PauliBitString(PauliString(pbs1))
    @test ps1 == PauliString(PauliBitString(ps1))

    ps1 = PauliString("YZZXYIIZXXYIZYIZZYYYIXYIYZZXYIIZIXYIYZZXYIIZZIXYIYZZXYIIZXXYIZYIZZYIXYIIXYIYZZXYIIZXXYIZYIZZYIXYI")
    ps2 = PauliString("XYYZYYZIYYYZYYZYYZZZYYZYYYYZYYZYYYZYYYYZYYZYYYYZYYYYZYYZYYZYYXIZYZZYYZZYYZYYYYZYYZYYZYYXIZYZZYYZZ")

    pbs1 = PauliBitString(ps1, T=UInt128)
    pbs2 = PauliBitString(ps2, T=UInt128)
    #@btime UnitaryPruning.commute($ps1, $ps2)
    #@btime UnitaryPruning.commute($pbs1, $pbs2)
    @test UnitaryPruning.commute(ps1, ps2) == UnitaryPruning.commute(pbs1, pbs2)

    results = []
    for i in 1:30
        for j in 1:30
            pbs1 = PauliBitString{UInt,20}(i,j)
            pbs2 = PauliBitString{UInt,20}(2*j,3*i)
            ps1 = PauliString(pbs1)
            ps2 = PauliString(pbs2)

            #println(ps1)
            #println(ps2)
            #println(UnitaryPruning.commute(ps1,ps2))
            #println(pbs1)
            #println(pbs2)
            #println(UnitaryPruning.commute(pbs1,pbs2))
            push!(results, UnitaryPruning.commute(pbs1,pbs2) == UnitaryPruning.commute(ps1,ps2))
        end
    end
    @test all(results)

    #println(" warntype function")
    #@code_warntype UnitaryPruning.commute(pbs1, pbs2) 
    #println(" warntype inline")
    #@code_warntype iseven(count_ones(pbs1.x & pbs2.z) - count_ones(pbs1.z & pbs2.x))
    #@btime iseven(count_ones($pbs1.x & $pbs2.z) - count_ones($pbs1.z & $pbs2.x))
    #@btime UnitaryPruning.commute($pbs1, $pbs2) 

    i = 5.6
    j = 982.1
    @btime $i<$j
end
