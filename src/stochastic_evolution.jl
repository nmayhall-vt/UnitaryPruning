using Random

"""
    stochastic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; nsamples=1000) where {N, P<:Pauli}

Stochastically sample the leaves of the binary tree where the probability of observing a leaf is proportional to its coefficient
"""
function stochastic_pauli_rotations(generators::Vector{Pauli{N}}, angles, o::Pauli{N}, ket; nsamples=1000) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)

    # 
    # Precompute the trig values to avoid doing so within the loop over samples
    scales = sin.(angles) .+ cos.(angles) 
    bias = sin.(angles) ./ scales
       
    # @btime stochastic_pauli_rotations_walk($generators, $o, $ket, $bias, $scales)

    #
    # collect our results here...
    expval = zeros(ComplexF64, nsamples)
  
    # population = PauliSum(N)
    population = Dict{Pauli{N},UInt}()
    
    # 
    # loop over number of samples for this run 
    for s in 1:nsamples
  
        #
        # initialize data for new walk down tree
        scale = 1.0
        nt = length(generators)
        oi = o

        # 
        # Loop through the generators in reverse, i.e., from the operator to the state
        for t in 1:nt
            g = generators[t]
           
            # @printf(" %12.8f %s\n", bias[t], string(g))
            #
            # First check to see if the current generator commutes with our current operator
            commute(oi.pauli, g.pauli) == false || continue

            #
            # Our bias goes from 0 -> 1 depending on if we should branch toward cos or sin respectively
            branch = rand() < bias[t]

            # if branch is true, we consider sin branch, else consider cos
            if branch
                # sin branch
                oi = g * oi    # multiply the pauli's
                oi = Pauli{N}((oi.θ + 1)%4, oi.pauli)             # multiply the sin branch pauli by 1im
            end
            scale *= scales[t]          # update the path scale with the current cos(a)+sin(a)
        end
        # sum!(population, oi)
        # population[oi] = get!(population, oi, 1)
        if haskey(population, oi)
            population[oi] += 1
        else
            population[oi] = 1
        end 
        expval[s] = scale * PauliOperators.expectation_value(oi, ket)
    end
    
    return expval, population
end




