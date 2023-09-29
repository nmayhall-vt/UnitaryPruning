using LinearAlgebra
using PushVectors
using Statistics 

"""
loop over each cos branch
"""
function spd(generators::Vector{Pauli{N}}, angles, o::Pauli{N}, ket ; thresh=1e-3, max_depth   =   100) where {T,N,P<:Pauli}

    nt = length(generators)

    vcos = cos.(angles)
    vsin = sin.(angles)

    # display(vcos)
    # display(vsin)
    path = zeros(Int8, nt+1)


    # oi, hi, ansatz_layer, branch_direction
    stack = Stack{Tuple{typeof(o), Float64, Int, Int8}}()  
    #stack = PushVector{Tuple{typeof(o), Float64, Int, Int, Int8}}()  
   
    err_estimate = 0.0

    # (operator, scale, layer, depth, branch_direction)
    push!(stack, (o, 1.0, 1, 0)) 

    final_value::Vector{T} = [0.0]
    final_paths::Vector{Int} = [0, 0]
    
    thresh2 = thresh*thresh

    iter = 0
    # stack contains a list of branch points from which we want to explore all left (cos) branches
    # each time we create a branch we add the right (sin) root to the stack
    while length(stack) > 0


        o, scale, layer, branch_direction = pop!(stack)
      
        
        layer_curr = layer

        for layer_idx in layer:nt
            layer_curr = layer_idx

            path[layer_idx] = branch_direction

            g = generators[layer_idx]
            if commute(g, o)
                path[layer_idx+1] = Int8(0)
                continue
            end


            # @printf(" %12.8f %12.8f\n", scale_l, scale_r)

            # should we branch right?
            # if abs(scale_r) > thresh && depth <= max_depth
            scale_r = scale * vsin[layer_idx]
            if abs(scale_r) > thresh
                or = multiply(o, g)
                or.θ = (or.θ + 1) % 4
                push!(stack, (or, scale_r, layer_idx + 1, -1))
            else
                err_estimate += abs(scale_r)
            end

            # left branch
            scale = scale * vcos[layer_idx]
            path[layer_idx+1] = 1
        
        end
            
        if layer_curr == nt 
            _found_leaf(ref_state, final_value, o, scale, path, final_paths)
            continue
        end


        # @printf(" gidx %4i  %12.8f %12.8f %s %s\n", layer_idx, scale_l, scale_r, string(o), string(g))
    end

    return final_value[1], final_paths[1], final_paths[2], err_estimate 
end



"""
loop over each cos branch
"""
function spd(ref_state::Vector{Bool}, o::PauliBoolVec{N}, generators::Vector{P}, angles::Vector{T};
    thresh      =   1e-12,
    max_depth   =   4) where {T,N,P<:Pauli}

    nt = length(generators)

    vcos = cos.(angles)
    vsin = sin.(angles)

    # display(vcos)
    # display(vsin)
    path = zeros(Int8, nt+1)


    # oi, hi, ansatz_layer, branch_direction
    stack = Stack{Tuple{typeof(o), Float64, Int, Int, Int8}}()  
    #stack = PushVector{Tuple{typeof(o), Float64, Int, Int, Int8}}()  
   
    err_estimate = 0.0

    # (operator, scale, layer, depth, branch_direction)
    push!(stack, (o, 1.0, 1, 0, 0)) 

    final_value::Vector{T} = [0.0]
    final_paths::Vector{Int} = [0, 0]
    
    thresh2 = thresh*thresh

    iter = 0
    # stack contains a list of branch points from which we want to explore all left (cos) branches
    # each time we create a branch we add the right (sin) root to the stack
    while length(stack) > 0


        # iter += 1
        # iter < 1000 || break
        o, scale, layer_idx, depth, branch_direction = pop!(stack)
        
        if layer_idx == nt+1 
            _found_leaf(ref_state, final_value, o, scale, path, final_paths)
            continue
        end
       
        
        path[layer_idx] = branch_direction

        g = generators[layer_idx]
        if commute(g,o)
            path[layer_idx+1] = Int8(0) 
            push!(stack, (o, scale, layer_idx + 1, depth, 1))
            continue
        end
            
        scale_l = scale * vcos[layer_idx]
        scale_r = scale * vsin[layer_idx]

        # @printf(" %12.8f %12.8f\n", scale_l, scale_r)
        
        # should we branch right?
        if abs(scale_r) > thresh && depth <= max_depth
            or = multiply(o, g)
            or.θ = (or.θ + 1) % 4
            push!(stack, (or, scale_r, layer_idx + 1, depth + 1, -1))
        else
            err_estimate += abs(scale_r)
        end
        
        # should we branch left?
        if abs(scale_l) > thresh
            push!(stack, (o, scale_l, layer_idx + 1, depth, 1))
        else
            err_estimate += abs(scale_l)
        end


        # @printf(" gidx %4i  %12.8f %12.8f %s %s\n", layer_idx, scale_l, scale_r, string(o), string(g))
    end

    return final_value[1], final_paths[1], final_paths[2], err_estimate 
end



function _found_leaf(ref_state, energy::Vector{T}, o::PauliBoolVec{N}, scale, path::Vector{Int8}, npaths::Vector{Int}) where {T,N}
    if is_diagonal(o)

        energy[1] += scale * expectation_value_sign(o, ref_state)
        npaths[1] += 1
    else
        npaths[2] += 1
    end
end