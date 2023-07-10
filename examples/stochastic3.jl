using UnitaryPruning
using Plots
using Statistics
using Printf
using Random
using LinearAlgebra

function run(;α=.01, k=10, nsamples=10000, seed=1)
    
    N = 20  # number of qubits

    generators = Vector{PauliBoolVec{N}}([])
    parameters = Vector{Float64}([])

    # Loop over trotter steps
    for ki in 1:k
        
        ## ZZ layer
        # e^{i π/2 P2} e^{i π P1 /2}|ψ>
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
    o = PauliBoolVec(N, Z=[1,2,3,4])


    ket = zeros(Bool, N)
    Random.seed!(seed)
    e = zeros(ComplexF64, nsamples)
    # e2 = zeros(ComplexF64, nsamples)
    
    ei, _ = UnitaryPruning.stochastic_pauli_dynamics_run(generators, parameters, o, ket)
    e[1] = ei 
    
    for i in 2:nsamples
        ei, _ = UnitaryPruning.stochastic_pauli_dynamics_run(generators, parameters, o, ket)
        e[i] = (e[i-1]*(i-1) + ei)/i
        # e2 += ei*ei
    end
    # @printf(" Mean: %12.8f Stdev: %12.8f\n", mean(results), std(results))

    # out = [results[1]]
    # for (idx,i) in enumerate(results)
    #     idx > 1 || continue
    #     push!(out, (out[idx-1]*(idx-1)+i)/idx)
    # end
    # e .= e ./ nsamples

    # e2 = e2 / nsamples
    # return e, e2 - e*e 
    return e
end

# e = run(;α=.01, k=100, nsamples=10000)

nruns = 100
nsamples = 1000

e_avgs = [];
e_vals = [];
v_vals = [];
for i in 0:16
    e = []
    etraj = []
    for runi in 1:nruns
        ei = run(α=i*π/32, k=5, nsamples=nsamples, seed=runi)
        push!(e, real(ei[end]))
        # push!(e, mean(real(ei)))
        # push!(v_vals, v)
        push!(etraj, ei)
    end
    push!(e_vals, etraj)
    avg = mean(e)
    push!(e_avgs, mean(etraj) )
    std = sqrt(sum(e.*e)) - avg*avg
    @printf(" avg: %12.8f var: %12.8f\n", avg, std)
    #@printf(" avg: %12.8f %12.8fi var: %12.8f %12.8fi\n", real(avg), imag(avg), real(std), imag(std))
    push!(v_vals, std)

    plot(real(e_vals[i+1]), ylim=[-1,1], legend=false, color=:grey, alpha=.1)
    plot!(real(e_avgs[i+1]), color=:black)
    savefig(@sprintf "%05i.pdf" i)
end


# plot(vals, marker = :circle)