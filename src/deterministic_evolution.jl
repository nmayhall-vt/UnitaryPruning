using Random
using BenchmarkTools
"""
    deterministic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; nsamples=1000) where {N, P<:Pauli}


"""
function deterministic_pauli_rotations(generators::Vector{Pauli{N}}, angles, o::Pauli{N}, ket ; thres=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    
    vcos = cos.(angles)
    vsin = sin.(angles)

    # collect our results here...
    expval = zero(ComplexF64)


    opers = PauliSum{N}(Dict(o.pauli=>1.0*(1im)^o.θ))#Dict{Tuple{Int128, Int128}, Complex{Float64}}((o.pauli.z,o.pauli.x)=>1.0*(1im)^o.θ)
        
    
    for t in 1:nt

        g = generators[t]
        branch_opers = PauliSum(N)#Dict{Tuple{Int128, Int128}, Complex{Float64}}()

#        sizehint!(branch_opers, 1000)
        for (key,value) in opers.ops
            
            oi = key
            if commute(oi, g.pauli)
                if haskey(branch_opers, oi)
                    branch_opers[oi] += value
                else
                    branch_opers[oi] = value
                end
                continue
            end
            if abs(value) > thres #If greater than threshold then split the branches
                # cos branch
                coeff = vcos[t] * value
                           
                if haskey(branch_opers, oi) # Add operator to dictionary if the key doesn't exist
                    branch_opers[oi] += coeff
                else
                    branch_opers[oi] = coeff #Modify the coeff if the key exists already
                end

                # sin branch
                coeff = vsin[t] * 1im * value
                oi = Pauli{N}(0, oi)
                oi = g * oi    # multiply the pauli's

                if haskey(branch_opers, oi.pauli) # Add operator to dictionary if the key doesn't exist
                    branch_opers[oi.pauli] += coeff * (1im)^oi.θ
                else
                    branch_opers[oi.pauli] = coeff * (1im)^oi.θ #Modify the coeff if the key exists already
                end
            end
        end
        opers = deepcopy(branch_opers) # Change the list of operators to the next row

    end

    for (key,value) in opers.ops
        oper = Pauli(UInt8(0), key)
        expval += value*PauliOperators.expectation_value_sign(oper, ket)
    end
    
    return expval
end




