using Distributed
@everywhere begin
    using UnitaryPruning
    using Plots
    using Statistics
    using Printf
    using Random
    using LinearAlgebra
    using SharedArrays
end

@everywhere function get_unitary_sequence_1D(o::PauliBitString{N}; α=.01, k=10) where N
    generators = Vector{PauliBitString{N}}([])
    parameters = Vector{Float64}([])
    # print("alpha", α, "\n")
    # Loop over trotter steps
    for ki in 1:k
        ## ZZ layer
        # e^{i π/2 P2} e^{i π P1 /2}|ψ>
        for i in 1:N-1
            pi = PauliBitString(N, Z=[i, i + 1])
            push!(generators, pi)
            push!(parameters, π/2)
        end
        #pbc 
        pi = PauliBitString(N, Z=[N, 1])
        push!(generators, pi)
        push!(parameters, π/2)

        ## X layer
        # e^{i αn (-X) / 2}
        for i in 1:N
            pi = PauliBitString(N, X=[i])
            pi += 2 # this accounts for the fact that the papers have -X and positive ZZ
            push!(generators, pi)
            push!(parameters, α)
        end
    end

    return generators, parameters
end

# @everywhere function get_unitary_sequence_1D(o::PauliBoolVec{N}; α=.01, k=10) where N
   
#     generators = Vector{PauliBoolVec{N}}([])
#     parameters = Vector{Float64}([])

#     # Loop over trotter steps
#     for ki in 1:k
        
#         ## ZZ layer
#         # e^{i π/4 P2} e^{i π P1 /4}|ψ>
#         for i in 1:N-1
#             pi = PauliBoolVec(N, Z=[i, i + 1])
#             push!(generators, pi)
#             push!(parameters, π/2)
#         end
            
#         #pbc 
#         pi = PauliBoolVec(N, Z=[N, 1])
#         push!(generators, pi)
#         push!(parameters, π/2 )

#         ## X layer
#         # e^{i αn Pn / 2}
#         for i in 1:N
#             pi = PauliBoolVec(N, X=[i])
#             pi += 2 # this accounts for the fact that the papers have -X and positive ZZ
#             push!(generators, pi)
#             push!(parameters, α)
#         end
#     end

#     return generators, parameters
# end

@everywhere function compute_run(generators, parameters, o, ket, N)
#    Random.seed!(seed)
#    rolling_avg = zeros(nsamples)

    samples = UnitaryPruning.deterministic_pauli_rotations(generators, parameters, o, ket, N)

    samples .= real.(samples)

    rolling_avg[1] = samples[1]
    
    for i in 2:nsamples
        rolling_avg[i] = (rolling_avg[i-1]*(i-1) + samples[i])/i
    end
    return rolling_avg, rolling_avg.*rolling_avg 
end


function run(; N=6)
   
    #
    # Uncomment to use bitstrings
    #
    ket = BasisState(N, 0) 
    o = PauliBitString(N, Z=[5])
    o = PauliBitString(N, Y=[1])
    # o = PauliBitString(N, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33])
    o = PauliBitString(N, Z=[1])

    #
    # Uncomment to use boolvecs
    #
    # ket = zeros(Bool,N)
    # o = PauliBoolVec(N, Z=[1, 2, 3, 4, 5, 6])
    # o = PauliBoolVec(N, Z=[1, 2])
    # o = PauliBoolVec(N, Y=[1])
    # o = PauliBoolVec(N, Z=[1])

    final_vals_stoc = []
    final_vals_errs = []

    for i in [(i-1)*2 for i in 1:16]

        
        #
        # Uncomment the following to do a serial run
        #
        generators, parameters = get_unitary_sequence_1D(o, α=i * π / 32, k=2)
        # for g in generators
        #     print(g)
        # end
        # a, b = compute_run(generators, parameters, o, ket, N)
    
        e = UnitaryPruning.deterministic_pauli_rotations(generators, parameters, o, ket, N)


        @printf(" α: %4i e: %12.8f+%12.8fi\n", i, real(e), imag(e))

    end

    # plot(final_vals_stoc, ribbon=final_vals_errs, marker=:circle)
    # plot!(final_vals)
    # savefig("plot.pdf")
    return final_vals_stoc, final_vals_errs
end

@time v,e = run(N=6)
