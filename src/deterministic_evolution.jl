using Random

"""
    deterministic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; nsamples=1000) where {N, P<:Pauli}


"""
function deterministic_pauli_rotations(generators::Vector{P}, angles, o::P, ket) where {N, P<:Pauli}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    print("nt", nt)

    # collect our results here...
    expval = zeros(ComplexF64,1)


    # Loop through the generators in reverse, i.e., from the operator to the state
    branch_opers = Dict{Tuple{Int128, Int128}, Float64}()
    opers = Dict((o.z,o.x)=>1.0)

    for t in reverse(1:nt)
#        oi = o
        branch_opers=[]
        g = generators[t]
        for oi in opers
        # First check to see if the current generator commutes with our current operator

            commute(oi, g) == false || continue

            if coeff > thres  #If greater than threshold then split the branches
                # cos branch
                coeff = cos.(angles[t])
                if !haskey(branch_opers, oi) 
                    branch_opers[(oi.z, oi.x)] = coeff * (1im)^oi.θ
                else
                    branch_opers[(oi.z, oi.x)] += coeff * (1im)^oi.θ
                end

                # sin branch
                oi = multiply(g, oi)    # multiply the pauli's
                oi = oi + 1             # multiply the sin branch pauli by 1im
                if !haskey(branch_opers, oi) 
                    branch_opers[(oi.z, oi.x)] = coeff * (1im)^oi.θ
                else
                    branch_opers[(oi.z, oi.x)] += coeff * (1im)^oi.θ
                end
             end
         end
         opers = branch_opers # Change the list of operators to the next row
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
