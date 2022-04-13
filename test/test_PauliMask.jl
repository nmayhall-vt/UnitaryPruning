using UnitaryPruning
using Test
using BenchmarkTools 

@testset "PauliMask" begin
    ps1 = PauliString("ZIYIZIYIY")

    ps2 = PauliString("XIYZYIXXY")

    @test UnitaryPruning.commute(ps1,ps2) == false

    pm1 = PauliMask(ps1)
    pm2 = PauliMask(ps2)

    @test is_diagonal(pm1) == false
    @test is_diagonal(PauliMask(PauliString("IZZIZIZZI"))) == true
    
    println(commute(pm1, pm2)) 
    println(commute(ps1, ps2))
    @test commute(pm1, pm2) == commute(ps1, ps2)
    
    @test pm1 == PauliMask(ps1)
    @test pm2 == PauliMask(ps2)

    @test ps1 == PauliString(pm1)
    @test ps2 == PauliString(pm2)

    println(string(ps1))
    println(string(ps2))
    println(commutator(ps1,ps2))
    println(commutator(pm1,pm2))
    #@btime commute($ps1,$ps2)
    #@btime commute($pm1,$pm2)
    #@btime commutator($ps1,$ps2)
    #@btime commutator($pm1,$pm2)
    #@code_warntype commutator(pm1,pm2)
    
    results = []
    for i in 1:3
        for j in 1:3
            pbs1 = PauliBitString{UInt,20}(i,j)
            pbs2 = PauliBitString{UInt,20}(2*j,3*i)
            ps1 = PauliString(pbs1)
            ps2 = PauliString(pbs2)

            pm1 = PauliMask(ps1)
            pm2 = PauliMask(ps2)
            
            push!(results, commute(pm1,pm2) == commute(ps1,ps2))
            
            if commute(ps1, ps2) == false
                phase1a, pm12 = commutator(pm1,pm2)
                phase1b, ps12 = commutator(ps1,ps2)
          
                @test  phase1a == phase1b
                @test  ps12 == PauliString(pm12)
                push!(results, phase1a == phase1b)
            end
        end
    end
    @test all(results)
end

