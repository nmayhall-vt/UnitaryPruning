using UnitaryPruning
using Plots
using Statistics
using Printf
using Random
using LinearAlgebra


function get_processor(N=127, k=10, α=.01, sequences=[[1,13]],#[19,32]],#,[38,51],[57,70],[76,89],[95,108],[114,126]],
                       bridges=[])#[[1,15,19],[5,16,23],[9,17,27],[13,18,31]])#,[33,37,52],[29,36,48],[25,35,44],[21,34,40],
#                                [38,53,57],[42,54,61],[46,55,65],[50,56,69],[71,75,90],[67,74,86],[63,73,82],[59,72,78],
#                                [76,91,95],[80,92,99],[84,93,103],[88,94,107],[97,110,115],[101,111,119],[105,112,123],
#                                [109,113,127]])

    need_bridges=false

    generators = Vector{PauliBoolVec{N}}([])
    parameters = Vector{Float64}([])

    # Loop over trotter steps
    for ki in 1:k        
        # e^{i π/2 P2} e^{i π P1 /2}|ψ>
        ## ZZ layer
        for qubit in sequences
            for i in qubit[1]:qubit[2]
                pi = PauliBoolVec(N, Z=[i, i + 1])
                push!(generators, pi)
                push!(parameters, π)
            end
        end
        #bridges
        if need_bridges
            for link in bridges
                pi = PauliBoolVec(N, Z=[link[1], link[2]])
                push!(generators, pi)
                push!(parameters, π )
                pi = PauliBoolVec(N, Z=[link[2], link[3]])
                push!(generators, pi)
                push!(parameters, π )
            end
        else #PBC
            pi = PauliBoolVec(N,Z=[1])
            push!(generators, pi)
            push!(parameters, π )
        end   
        ## X layer
        # e^{i αn Pn / 2}
        for i in 1:N
            pi = PauliBoolVec(N, X=[i])
            push!(generators, pi)
            push!(parameters, α)
        end
    end
    return(generators,parameters)
end

function run(;α=.01, k=10)
    N = 14

    generators,parameters = get_processor(N,k,α)
    # o = PauliBoolVec(N,X=[1],Y=[2],Z=[3])

    #Mz
    # o = PauliBoolVec(N, Y=[1], Z=[2,3,4])
    o = PauliBoolVec(N, Z=[1])
    o_mat = to_matrix(o)
    # for i in 2:N
    #     o .+= to_matrix(PauliBoolVec(N, Z=[i]))
    # end


    U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    m = diag(U'*o_mat*U)
    # println(" expectation values:")
    # display(m[1])
   
    return real(m[1])
end

# run(;α=.01, k=30)

final_vals = []
for i in 0:1
    ei = run(α=i*π/32, k=1)
    push!(final_vals, real(ei))
    @printf(" α: %4i val: %12.8f\n", i, ei)
end


# vals = [];
# for i in 1:40
#     push!(vals, run(α=0.1, k=i))
# end

plot(final_vals, marker = :circle)
