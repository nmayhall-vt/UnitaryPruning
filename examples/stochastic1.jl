using UnitaryPruning
using Plots
using Statistics
using Printf
using Random
using LinearAlgebra

function run(;α=.01, k=10)
    N = 6 
    generators = Vector{PauliBoolVec{N}}([])
    parameters = Vector{Float64}([])

    # Loop over trotter steps
    for ki in 1:k
        # e^{i π/2 P2} e^{i π P1 /2}|ψ>
        ## ZZ layer
        for i in 1:N-1
            pi = PauliBoolVec(N, Z=[i, i + 1])
            push!(generators, pi)
            push!(parameters, π)
        end
            
        #pbc 
        pi = PauliBoolVec(N, Z=[N, 1])
        push!(generators, pi)
        push!(parameters, π )

        ## X layer
        # e^{i αn Pn / 2}
        for i in 1:N
            pi = PauliBoolVec(N, X=[i])
            push!(generators, pi)
            push!(parameters, α)
        end
    end

    # o = PauliBoolVec(N,X=[1],Y=[2],Z=[3])

    #Mz
    o = PauliBoolVec(N, Z=[1])
    o_mat = to_matrix(o)
    # for i in 2:N
    #     o .+= to_matrix(PauliBoolVec(N, Z=[i]))
    # end


    U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    
    # Ui = UnitaryPruning.build_time_evolution_matrix([generators[2]], [parameters[2]])
    # Uj = exp(to_matrix(generators[2]) .* 1im .* parameters[2])

    # display(norm(Ui - Uj))

    # out = [mean(results[1:i]) for i in 1:length(results)]
    # plot(real(results))
    m = diag(U'*o_mat*U)
    println(" expectation values:")
    display(m[1])
   
    # return real(m[1])
    
    ket = zeros(Bool, N)
    Random.seed!(2)
    results = Vector{ComplexF64}([])
    for i in 1:1000
        e, _ = UnitaryPruning.stochastic_pauli_dynamics_run(generators, parameters, o, ket)
        # display(path)
        push!(results, e)
    end
    @printf(" Mean: %12.8f Stdev: %12.8f\n", mean(results), std(results))

    out = [results[1]]
    for (idx,i) in enumerate(results)
        idx > 1 || continue
        push!(out, (out[idx-1]*(idx-1)+i)/idx)
    end

    return real(m[1]), out 
end

run(;α=.01, k=30)


# vals = [];
# for i in 1:40
#     push!(vals, run(α=0.1, k=i))
# end

# plot(vals, marker = :circle)