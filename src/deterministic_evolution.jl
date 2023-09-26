using Random
using BenchmarkTools
"""
    deterministic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; threshold) where {N, P<:Pauli}


"""
function deterministic_pauli_rotations_BFS(generators::Vector{Pauli{N}}, angles, o::Pauli{N}, ket ; thres=1e-3) where {N}

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



function deterministic_pauli_rotations_DFS(generators::Vector{Pauli{N}}, angles, o::Pauli{N}, ket ; thres=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    
    vcos = cos.(angles)
    vsin = sin.(angles)

    # collect our results here...
    expval = zero(ComplexF64)


    opers = ScaledPauliVector{N}([ScaledPauli(o)])
        
    
    for t in 1:nt

        g = generators[t]
        branch_opers = ScaledPauliVector(N)

#        sizehint!(branch_opers, 1000)
        for sp in opers
            
            oi = sp.pauli
            if commute(oi, g.pauli)
                push!(branch_opers, sp)
                continue
            end
            if abs(sp.coeff) > thres #If greater than threshold then split the branches
                # cos branch
                coeff = vcos[t] * sp.coeff
                push!(branch_opers, ScaledPauli{N}(coeff, oi))
                           
                # sin branch
                coeff = vsin[t] * 1im * sp.coeff
           
                oi = Pauli{N}(0, oi)
                oi = g * oi    # multiply the pauli's
                coeff = coeff * (1im)^oi.θ
                push!(branch_opers, ScaledPauli{N}(coeff, oi.pauli))
            end
        end
        opers = deepcopy(branch_opers) # Change the list of operators to the next row

    end

    for sp in opers
        oper = Pauli(UInt8(0), sp.pauli)
        expval += sp.coeff*PauliOperators.expectation_value_sign(oper, ket)
    end
    
    return expval
end




