using Distributed 


testcommute(p1::FixedPhasePauli, p2::FixedPhasePauli) = iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 
testcommute(p1::Pauli, p2::Pauli) = iseven(count_ones(p1.pauli.x & p2.pauli.z) - count_ones(p1.pauli.z & p2.pauli.x)) 

"""
    iterate_dfs!(ket, o::P, h::T, generators::Vector{Vector{P}}, 
                      ansatz_par::Vector{T}; thresh=1e-12, max_depth=20) where {T,N, P<:Pauli}

"""
function dfs_evolution(generators::Vector{Pauli{N}}, angles, h, o::Pauli{N}, ket ; thresh=1e-3, max_depth=4) where {N}

    vcos = cos.(angles)
    vsin = sin.(angles)
  
    # Operator, scale
    stack = Stack{Tuple{Pauli{N},ComplexF64,Int}}()  

    n_gen = length(generators)

    for op_idx in 1:n_gen
        if commute(generators[op_idx], o) == false
            push!(stack, (o, h, op_idx))
            break
        end
    end 

    my_expval::Vector{ComplexF64} = [0.0]
    l2::Vector{Float64} = [0, 0]

    while length(stack) > 0
        oi, hi, op_idx = pop!(stack)
            
    
        # if op_idx == n_gen+1
        #     _found_leaf(ket, my_expval, oi, hi, l2)
        # elseif sqrt(abs2(hi)) > thresh
        
        if sqrt(abs2(hi)) > thresh
            for op_idx2 in op_idx:n_gen
                g = generators[op_idx2]
                
                if testcommute(g, oi) == false
            
                    # sin branch
                    or = g*oi
                    hr = 1im * hi * vsin[op_idx2]

                    if op_idx2 == n_gen
                        _found_leaf(ket, my_expval, or, hr, l2)
                    else
                        push!(stack, (or, hr, op_idx2+1))
                    end
            
                    # cos branch
                    hl = hi * vcos[op_idx2]
                    if op_idx2 == n_gen
                        _found_leaf(ket, my_expval, oi, hl, l2)
                    else
                        push!(stack, (oi, hl, op_idx2+1))
                    end

                    break
                end 
                if op_idx2 == n_gen
                    _found_leaf(ket, my_expval, oi, hi, l2)
                end
            end 
                
        end
    end
    return my_expval[1], l2[1], l2[2]
    # return my_expval[1], l2[1], l2[2] 
end


function dfs_evolution2(generators::Vector{Pauli{N}}, angles, h, o::Pauli{N}, ket ; thresh=1e-3, max_depth=4) where {N}

    vcos = cos.(angles)
    vsin = sin.(angles)
    depth = 0
  
    # Operator, scale, Layer, depth
    stack = Stack{Tuple{Pauli{N},ComplexF64,Int,Int}}()  

    n_gen = length(generators)

    push!(stack, (o,h,1,1)) 

    my_expval::Vector{ComplexF64} = [0.0]
    l2::Vector{Float64} = [0, 0]

    while length(stack) > 0
        oi, hi, op_idx, depth = pop!(stack)
    
        if op_idx == n_gen+1
            _found_leaf(ket, my_expval, oi, hi, l2 )
        elseif abs(hi) < thresh
            # _found_leaf(ket, my_expval, oi, hi, l2 )
        else
            g = generators[op_idx]
            
            if commute(g,oi)
                push!(stack, (oi, hi, op_idx+1, depth))
            else
                # if depth >= max_depth
                #     _found_leaf(ket, my_expval, oi, hi, l2 )
                # end

                # right branch
                or = g*oi
                hr = 1im * hi * vsin[op_idx]

                push!(stack, (or, hr, op_idx+1, depth+1))

                # left branch
                hl = hi * vcos[op_idx]
                push!(stack, (oi, hl, op_idx+1, depth))
            end
        end
    end
    return my_expval[1], l2[1], l2[2]
    # return my_expval[1], l2[1], l2[2] 
end


function _found_leaf(ket, energy::Vector, o, h, l2)
    # error("nick")
    if is_diagonal(o)
        sign = expectation_value(o, ket) 

        energy[1] += sign*h
        l2[1] += h*h' 
    else
        l2[2] += h*h' 
    end
end