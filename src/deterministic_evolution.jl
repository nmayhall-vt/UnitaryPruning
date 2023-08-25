using Random
using BenchmarkTools
"""
    deterministic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; nsamples=1000) where {N, P<:Pauli}


"""
function deterministic_pauli_rotations(generators::Vector{PauliBitString{N}}, angles, o::PauliBitString{N}, ket ; thres=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    
    vcos = cos.(angles)
    vsin = sin.(angles)

    # collect our results here...
    expval = zero(ComplexF64)


    # Loop through the generators in reverse, i.e., from the operator to the state
    opers = Dict{Tuple{Int128, Int128}, Complex{Float64}}((o.z,o.x)=>1.0*(1im)^o.θ)
        
    
    for t in reverse(1:nt)

        g = generators[t]
        branch_opers = Dict{Tuple{Int128, Int128}, Complex{Float64}}()

        sizehint!(branch_opers, 1000)
        for (key,value) in opers
            
            oi = PauliBitString{N}(0,key[1],key[2])
            if commute(oi, g)
                if haskey(branch_opers, (oi.z,oi.x))
                    branch_opers[(oi.z, oi.x)] += (1im)^oi.θ * value
                else
                    branch_opers[(oi.z, oi.x)] = (1im)^oi.θ * value
                end
                continue
            end

            if abs(value) > thres #If greater than threshold then split the branches
                coeff = vcos[t] * value
            
               
                # cos branch
                if !haskey(branch_opers, (oi.z,oi.x)) # Add operator to dictionary if the key doesn't exist
                    branch_opers[(oi.z, oi.x)] = coeff * (1im)^oi.θ
                else
                    branch_opers[(oi.z, oi.x)] += coeff * (1im)^oi.θ #Modify the coeff if the key exists already
                end

                # sin branch
                coeff = vsin[t] * 1im * value
                oi = multiply(g, oi)    # multiply the pauli's

                if !haskey(branch_opers, (oi.z, oi.x)) # Add operator to dictionary if the key doesn't exist
                    branch_opers[(oi.z, oi.x)] = coeff * (1im)^oi.θ
                else
                    branch_opers[(oi.z, oi.x)] += coeff * (1im)^oi.θ #Modify the coeff if the key exists already
                end
            end
        end
        
        opers = deepcopy(branch_opers) # Change the list of operators to the next row

    end
    for (key,value) in opers
        oper = PauliBitString{N}(0, key[1], key[2])
        expval += value*expectation_value_sign(oper, ket)
    end
    
    return expval
end




