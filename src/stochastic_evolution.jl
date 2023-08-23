using Random

"""
    stochastic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; nsamples=1000) where {N, P<:Pauli}

Stochastically sample the leaves of the binary tree where the probability of observing a leaf is proportional to its coefficient
"""
function stochastic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; nsamples=1000) where {N, P<:Pauli}

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
        for t in reverse(1:nt)
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
                oi = multiply(g, oi)    # multiply the pauli's
                oi = oi + 1             # multiply the sin branch pauli by 1im
            end
            scale *= scales[t]          # update the path scale with the current cos(a)+sin(a)
        end
    
        expval[s] = scale * expectation_value_sign(oi, ket)
    end
    
    return expval
end





""""
    stochastic_pauli_rotations_run(generators::Vector{PauliBoolVec{N}}, angles, o::PauliBoolVec{N}) where N

Not yet used for anything
"""
function get_random_leaf(generators::Vector{PauliBoolVec{N}}, angles, o_in::PauliBoolVec{N}, ket) where N

    o = deepcopy(o_in)
    #
    # for a single pauli Unitary, Un = exp(-i θn Pn/2)
    # U' O U = cos(θ/2) O - i sin(θ/2) PO
    nt = length(generators)
    length(angles) == nt || throw(DimensionMismatch)

    scale = 1.0
    path = 1
    sin_angles = sin.(angles)
    cos_angles = cos.(angles)
  
    for t in reverse(1:nt)
        g = generators[t]
        commute(o,g) == false || continue

       

        if rand() < .5 
            # sin branch

            # rand() < sin_angles[t] || return 0.0, path
            # rand() < sin_angles[t] || continue 

            # println("sin: ", bias[t])
            # error("here")
            o = multiply(g,o)
            o.θ = (o.θ + 3) % 4

            scale *= sin_angles[t]

        else
            scale *= cos_angles[t]
        end
    end

    # return expectation_value_sign(o,ket), path
    return scale * expectation_value_sign(o,ket), path

end