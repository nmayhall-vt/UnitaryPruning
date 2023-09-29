using UnitaryPruning
using Test
using BenchmarkTools 
using PauliOperators
using LinearAlgebra

@testset "deterministic_pauli_rotations_BFS" begin

    N = 6
    ket = KetBitString(N, 0)
    
    op = [Pauli(N, X=[2,3], Y=[4], Z=[1,5]), Pauli(N, X=[3,4], Y=[1,2], Z=[5]), Pauli(N, X=[1], Y=[3,4], Z=[2,6]), Pauli(N, X=[1,2], Y=[6], Z=[3,4])]

    for o in op
        for i in 1:2:16
            α = i * π/32
            generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=2)
            e = UnitaryPruning.deterministic_pauli_rotations_BFS(generators, parameters, o, ket, thres=1e-4)
            o_mat = Matrix(o)
            U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
            m = diag(U'*o_mat*U)

            println(m[1], e)
            @test(abs(real(e)-real(m[1])) <= 1e-3)
        end
     end

end
