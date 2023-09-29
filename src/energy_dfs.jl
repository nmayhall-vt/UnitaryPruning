using Distributed 



"""
    iterate_dfs!(ket, o::P, h::T, generators::Vector{Vector{P}}, 
                      ansatz_par::Vector{T}; thresh=1e-12, max_depth=20) where {T,N, P<:Pauli}

"""
function dfs_evolution(generators::Vector{Pauli{N}}, angles, h, o::Pauli{N}, ket ; thresh=1e-3, max_depth=4) where {N}

    vcos = cos.(angles)
    vsin = sin.(angles)
    depth = 0
  
    # Operator, scale, Layer, depth
    stack = Stack{Tuple{Pauli{N},ComplexF64,Int,Int}}()  

    n_gen = length(generators)

    push!(stack, (o,h,1,1)) 

    my_expval::Vector{ComplexF64} = [0.0]
    my_paths::Vector{Int} = [0, 0]

    while length(stack) > 0
        oi, hi, op_idx, depth = pop!(stack)
    
        if op_idx == n_gen+1
            _found_leaf(ket, my_expval, oi, hi, my_paths )
        elseif abs(hi) < thresh
            _found_leaf(ket, my_expval, oi, hi, my_paths )
        else
            g = generators[op_idx]
            
            if commute(g,oi)
                push!(stack, (oi, hi, op_idx+1, depth))
            else
                # if depth >= max_depth
                #     _found_leaf(ket, my_expval, oi, hi, my_paths )
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
    return my_expval[1], my_paths[1], my_paths[2] 
end


function _found_leaf(ket, energy::Vector, o, h, paths::Vector{Int})
        if is_diagonal(o)
            sign = expectation_value(o, ket) 
    
            energy[1] += sign*h
            paths[1] += 1
        else
            paths[2] += 1
        end
    end