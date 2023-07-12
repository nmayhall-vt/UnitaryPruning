using Distributed
#@everywhere begin
using UnitaryPruning
using Plots
using Statistics
using Printf
using Random
using LinearAlgebra
using SharedArrays
#end

function get_unitary_sequence_1D(;α=.01, k=10, N=6)
    
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

    return generators, parameters
end


function run()
    
    N = 6

    # Operator
    o = PauliBoolVec(N, Z=[1])
    o = PauliBoolVec(N, Z=[1, 2])
    o = PauliBoolVec(N, Y=[1, 2, 3, 4, 5, 6])

    # State
    ref_state = zeros(Bool, N)

    final_vals_spd = []
    final_vals_err = []

    # for i in 8:8
    for i in 0:16
        generators, angles = get_unitary_sequence_1D(α=i * π / 32, k=5, N=N)
        
        # display(angles)
        # error("nick")
       
        # value, npath_nonzero, npath_zero, err = UnitaryPruning.spd2(ref_state, o, generators, angles, thresh=1e-4)
        value, npath_nonzero, npath_zero, err = UnitaryPruning.spd(ref_state, o, generators, angles, thresh=1e-4, max_depth=10)

        push!(final_vals_spd, value)
        push!(final_vals_err, err)
        @printf(" SPD = %12.8f Error < %12.8f # Diagonal Paths %12i  # Nondiagonal Paths %12i\n", value, err, npath_nonzero, npath_zero)
    end
    
    plot(final_vals_spd, ribbon=final_vals_err, marker=:circle)
    plot!(final_vals)
    savefig("plot.pdf")
end


run()