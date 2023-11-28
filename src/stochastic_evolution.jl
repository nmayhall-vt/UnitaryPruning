using Random

"""
    stochastic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; sample_size=1000) where {N, P<:Pauli}

Stochastically sample the leaves of the binary tree where the probability of observing a leaf is proportional to its coefficient
"""
function stochastic_pauli_rotations(generators::Vector{Pauli{N}}, angles, o::Pauli{N}, ket; sample_size=1000) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    
        nt = length(generators)
    length(angles) == nt || throw(DimensionMismatch)

    # 
    # Precompute the trig values to avoid doing so within the loop over samples
    bias = sin.(angles) .* sin.(angles) 
      

    #
    # sample = PauliSum(N)
    sample = Dict{Pauli{N},UInt}()
    
    # 
    # loop over number of samples for this run 
    for s in 1:sample_size
  
        #
        # initialize data for new walk down tree
        oi = deepcopy(o)

        # 
        # Loop through the generators in reverse, i.e., from the operator to the state
        for t in 1:nt
            g = generators[t]
           
            # @printf(" %12.8f %s\n", bias[t], string(g))
            #
            # First check to see if the current generator commutes with our current operator
            commute(oi, g) == false || continue

            #
            # Our bias goes from 0 -> 1 depending on if we should branch toward cos or sin respectively
            branch = rand() < bias[t]

            # if branch is true, we consider sin branch, else consider cos
            if branch
                # sin branch
                oi = oi * g     # multiply the pauli's
                oi = Pauli{N}((oi.θ + 1)%4, oi.pauli)             # multiply the sin branch pauli by 1im
            end
        end
        if haskey(sample, oi)
            sample[oi] += 1
        else
            sample[oi] = 1
        end 
    end
  
    effective_sample_size = 0

    # At this point, sample contains a list of up to 4*4^N Pauli's, which is larger than the full set of 
    # 4^N Pauli's. The reason is that each Pauli contains a phase of either 1,-1,i,-i. 
    H = PauliSum(N)
    for (oi, count) in sample
        prob_i = count / sample_size
        prob_i > 0 || throw(ErrorException)
        coeff = sqrt(prob_i)
        sum!(H, coeff*oi)
    end

    
    l2 = 0.0
    for (oi, coeff) in H.ops
        l2 += abs2(coeff)
    end
    l2 = sqrt(l2)
    H = H * (1/l2)
   
    # @show length(sample)
    # @show length(H.ops)
    return sample, H
end




