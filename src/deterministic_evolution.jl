using Random

"""
    deterministic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; nsamples=1000) where {N, P<:Pauli}


"""
function deterministic_pauli_rotations(generators::Vector{P}, angles, o::P, ket, N) where {P<:Pauli}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    thres = 10^-3
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
#    print("nt", nt)

    # collect our results here...
    expval = zeros(ComplexF64,1)


    # Loop through the generators in reverse, i.e., from the operator to the state
    opers = Dict((o.z,o.x)=>1.0*(1im)^o.θ)
    
    for t in reverse(1:nt)

        g = generators[t]
        branch_opers = Dict{Tuple{Int128, Int128}, Complex{Float64}}()
        for (key,value) in opers
        # First check to see if the current generator commutes with our current operator
#            print("oi", PauliBitString(key[1],key[2],N))
            oi = PauliBitString(key[1],key[2],N)

            println("commutes  ", commute(oi, g))
            println("oi",oi)
            println("g",g)
            commutes = commute(oi, g)
            if commutes
                branch_opers[(oi.z, oi.x)] = value
            end
            
            commutes == false || continue
            if (abs(real(value)) > thres || abs(imag(value)) > thres)#If greater than threshold then split the branches
                println("Above Threshold")
                # cos branch                
                coeff = cos.(angles[t])
                
                if !haskey(branch_opers, (oi.z,oi.x)) # Add operator to dictionary if the key doesn't exist
                    branch_opers[(oi.z, oi.x)] = coeff * (1im)^oi.θ
                    println("cos-without key", oi)
                else
                    branch_opers[(oi.z, oi.x)] += coeff * (1im)^oi.θ #Modify the coeff if the key exists already
                    println("cos-with key", oi)
                end

                # sin branch
                coeff = sin.(angles[t])
                oi = multiply(g, oi)    # multiply the pauli's
                oi = oi + 1             # multiply the sin branch pauli by 1im

                if !haskey(branch_opers, (oi.z,oi.x)) # Add operator to dictionary if the key doesn't exist
                    branch_opers[(oi.z, oi.x)] = coeff * (1im)^oi.θ
                    println("sin-without key", oi)
                else
                    branch_opers[(oi.z, oi.x)] += coeff * (1im)^oi.θ #Modify the coeff if the key exists already
                    println("sin-with key", oi)
                end
            end
        end
        if isempty(branch_opers)
            opers = opers
        else
            opers = branch_opers # Change the list of operators to the next row
        end 
        println("opers", opers)

    end
    println("Final", opers)
    for (key,value) in opers
        oper = PauliBitString(key[1],key[2],N)
#        println("oper  ", oper, is_diagonal(oper))
#        println("expval  ", expectation_value_sign(oper, ket))
        expval[1] += expectation_value_sign(oper, ket)
    end
    print(expval)
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
