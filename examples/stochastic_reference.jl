using UnitaryPruning
using Plots
using Statistics
using Printf
using Random
using LinearAlgebra

#
#   in the experiment, the circuit is 
#
#   exp(i θ/2 (-X)) exp(i π/4 ZZ)
#
function run(;α=.01, k=10)
    N = 6 
    generators = Vector{PauliBoolVec{N}}([])
    parameters = Vector{Float64}([])

    # Loop over trotter steps
    for ki in 1:k
        
        ## ZZ layer
        # e^{i θ/2 P2} e^{i θ/2 P1}|ψ>
        for i in 1:N-1
            pi = PauliBoolVec(N, Z=[i, i + 1])
            push!(generators, pi)
            push!(parameters, π/2)
        end
            
        #pbc 
        pi = PauliBoolVec(N, Z=[N, 1])
        push!(generators, pi)
        push!(parameters, π/2 )

        ## X layer
        # e^{i θ/2 (-P2)}
        for i in 1:N
            pi = PauliBoolVec(N, X=[i])
            pi += 2 # this accounts for the fact that the papers have -X and positive ZZ
            push!(generators, pi)
            push!(parameters, α)
        end
    end

    # o = PauliBoolVec(N,X=[1],Y=[2],Z=[3])

    #Mz
    o = PauliBoolVec(N, Y=[1], Z=[2,3,4])
    o = PauliBoolVec(N, Z=[1,2])
    o = PauliBoolVec(N, X=[1,2], Y=[3], Z=[5,6])
    o = PauliBoolVec(N, Z=[1,2,3,4,5,6])
    o = PauliBoolVec(N, Z=[1,2])
    o = PauliBoolVec(N, Y=[1])
    o = PauliBoolVec(N, X=[1], Y=[2,3])
    o = PauliBitString(N, Z=[1])
    
    o_mat = to_matrix(o)


    U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    
    m = diag(U'*o_mat*U)
    # println(" expectation values:")
    # display(m[1])
   
    return real(m[1])
end

# run(;α=.01, k=30)

final_vals = []
for i in 0:16
    ei = run(α=i*π/32, k=5)
    push!(final_vals, real(ei))
    @printf(" α: %4i val: %12.8f\n", i, ei)
end


# vals = [];
# for i in 1:40
#     push!(vals, run(α=0.1, k=i))
# end

plot(final_vals, marker = :circle)
