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

function run(;α=.01, k=10, N=6)


    # o = Pauli(N,X=[1],Y=[2],Z=[3])

    #Mz
#    o = Pauli(N, Y=[1], Z=[2,3,4])
#    o = Pauli(N, Z=[1,2])
#    o = Pauli(N, X=[1,2], Y=[3], Z=[5,6])
#    o = Pauli(N, Z=[1,2,3,4,5,6])
#    o = Pauli(N, Z=[2,5,7,8], Y=[1,3,4])
#    o = Pauli(N, X=[1], Y=[2,3])
#    o = Pauli(N, Z=[1,2])
    o = Pauli(N, Z=[1])
    
    generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=2)

    o_mat = Matrix(o)

    U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    
    m = diag(U'*o_mat*U)
    # println(" expectation values:")
    # display(m[1])
   
    return real(m[1])
end

# run(;α=.01, k=30)

final_vals = []
for i in [(i-1)*2 for i in 1:16]
    ei = run(α=i*π/32, k=2, N=3)
    push!(final_vals, real(ei))
    @printf(" α: %4i val: %12.8f\n", i, ei)
end


# vals = [];
# for i in 1:40
#     push!(vals, run(α=0.1, k=i))
# end

plot(final_vals, marker = :circle)
savefig("plot.pdf")


