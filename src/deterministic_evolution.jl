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


    # opers = PauliSum{N}(Dict(o.pauli=>1.0*(1im)^o.θ))#Dict{Tuple{Int128, Int128}, Complex{Float64}}((o.pauli.z,o.pauli.x)=>1.0*(1im)^o.θ)
    opers = PauliSum(o)
  
    n_ops = zeros(Int,nt)
    
    for t in 1:nt

        g = generators[t]
        branch_opers = PauliSum(N)

        sizehint!(branch_opers.ops, 1000)
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
        n_ops[t] = length(branch_opers)
        opers = deepcopy(branch_opers) # Change the list of operators to the next row

    end

    for (key,value) in opers.ops
        oper = Pauli(UInt8(0), key)
        expval += value*PauliOperators.expectation_value(oper, ket)
    end
   
    return expval, n_ops
end

function weight(p::FixedPhasePauli) 
    return  count_ones(p.x | p.z)
end

"""
    my_clip!(ps::PauliSum{N}; thresh=1e-16) where {N}

Clip small coefficients, with bias against higher weight pauli's
"""
function my_clip!(ps::PauliSum{N}; thresh=1e-16) where {N}
    filter!(p->abs(p.second)/weight(p.first) > thresh, ps.ops)
    # filter!(p->abs(p.second)/sqrt(weight(p.first)) > thresh, ps.ops)
end

"""
    normalize!(ps)

Normalize a `PauliSum`, `ps`
"""
function normalize!(ps)
    l2 = 0.0
    for (oi,coeff) in ps.ops
        l2 += abs2(coeff)
    end
    l2 = sqrt(l2)
    for (oi,coeff) in ps.ops
        ps.ops[oi] /= l2
    end
    return l2
end

"""
    bfs_evolution(generators::Vector{Pauli{N}}, angles, o::Pauli{N}, ket ; thres=1e-3) where {N}


"""
function bfs_evolution(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    
    vcos = cos.(angles)
    vsin = sin.(angles)



    o_transformed = deepcopy(o)
  
    n_ops = zeros(Int,nt)
    
    for t in 1:nt

        g = generators[t]

        sin_branch = PauliSum(N)

        for (oi,coeff) in o_transformed.ops
           
            abs(coeff) > thresh || continue


            if commute(oi, g.pauli) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
                # if count_ones(oj.pauli.x | oj.pauli.z) < max_weight   
                #     sum!(sin_branch, oj * vsin[t] * coeff * 1im)
                # end

            end
        end
        sum!(o_transformed, sin_branch) 
        # clip!(o_transformed, thresh=thresh)
        my_clip!(o_transformed, thresh=thresh)
        # normalize!(o_transformed) #this makes things worse
        n_ops[t] = length(o_transformed)
    end
    
    l1 = 0.0
    l2 = 0.0
    l4 = 0.0
    expval = zero(ComplexF64)
    avgval = 0.0
    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
        # if real(coeff) < 0.0
        #     avgval += abs2(coeff)*expectation_value(oi,ket)
        # else
        #     avgval += abs2(coeff)*expectation_value(oi,ket)
        # end
        l1 += abs(coeff)
        l2 += abs2(coeff)
        l4 += abs2(abs2(coeff))
    end
   
    entropy = 0.0
    for (oi,coeff) in o_transformed.ops
        pi = abs2(coeff)/l2
        entropy -= pi * log(2,pi) 
    end
   
    return expval, n_ops, l1, sqrt(l2), l4^.25, entropy, o_transformed
end


