using UnitaryPruning
using Plots
using Statistics
using Printf
using Random
using LinearAlgebra
using PauliOperators
#
#   in the experiment, the circuit is 
#
#   exp(i θ/2 (-X)) exp(i π/4 ZZ)
#

function run(;k=5, N=6)

    o = Pauli(N, Z=[1])
    o_mat = Matrix(o)

    e = []
    angles = []
    for i in  0:16
    # for i in [(i-1)*2 for i in 1:9]
        α=i*π/32
        
        generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=k)
        U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
        ei = diag(U'*o_mat*U)[1]

        push!(e, real(ei))
        push!(angles, α)
        
        @printf(" i: %4i α: %12.8f val: %12.8f\n", i, α, ei)
    end

    plot(angles, real(e), marker = :circle)
    xlabel!("Angle")
    savefig("plot_1d_n6_reference.pdf")
    
    return 
end

run(k=5, N=6)



