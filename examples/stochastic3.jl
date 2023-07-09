using UnitaryPruning
using Plots
using Statistics
using Printf
using Random
using LinearAlgebra

function run(;α=.01, k=10, nsamples=10000)
    N = 40 
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
    # o = PauliBoolVec(N, X=[13,29,31], Y=[9,30], Z=[8,12,17,28,32])
    o = PauliBoolVec(N, Z=[1])


    ket = zeros(Bool, N)
    Random.seed!(2)
    results = Vector{ComplexF64}([])
    e = 0.0
    for i in 1:nsamples
        ei, _ = UnitaryPruning.stochastic_pauli_dynamics_run(generators, parameters, o, ket)
        e += ei
    end
    # @printf(" Mean: %12.8f Stdev: %12.8f\n", mean(results), std(results))

    # out = [results[1]]
    # for (idx,i) in enumerate(results)
    #     idx > 1 || continue
    #     push!(out, (out[idx-1]*(idx-1)+i)/idx)
    # end
    e = e / nsamples

    return e 
end

# e = run(;α=.01, k=100, nsamples=10000)


vals = [];
for i in 1:16
    push!(vals, run(;α=i*π/32, k=6, nsamples=10000))
    println(vals[end])
end

# plot(vals, marker = :circle)