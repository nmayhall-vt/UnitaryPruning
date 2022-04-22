using Distributed 


"""
    compute_expectation_value_iter_parallel(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par; thresh=1e-8, max_depth=20)
"""
function compute_expectation_value_iter_parallel(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par; thresh=1e-8, max_depth=20)
#={{{=#
   
    energy, npath_nonzero, npath_zero = @distributed (.+) for hi in 1:length(ham_ops)
        #ei = expectation_value_sign(ham_ops[hi], ref_state) * ham_par[hi] 
        #e_hf += ei
    
        #energy::Vector{Float64} = [0.0]
        #if is_diagonal(ham_ops[hi])
        #    @printf("%20s %12.8f %12.8f\n", string(ham_ops[hi]), ei, ham_par[hi])
        #end
        #iterate_dfs!(ref_state, [0.0], [0,0], ham_ops[hi], ham_par[hi], ansatz_ops, ansatz_par, thresh=thresh, max_depth=max_depth)
        iterate_dfs!(ref_state, 
                         ham_ops[hi], ham_par[hi], 
                         ansatz_ops, ansatz_par, 
                         thresh=thresh, max_depth=max_depth)
    end
    #@printf(" E(HF)     = %12.8f\n", e_hf)
    #@printf(" E(cADAPT) = %12.8f\n", energy[1])
    #@printf(" %% Contributing Branches %12.4f %%  Tot: %i\n", paths[1]/sum(paths)*100, sum(paths[1]) )
    @printf(" E(cADAPT) = %12.8f  #Diagonal Paths %12i  #Nondiagonal Paths %12i\n", energy, npath_nonzero, npath_zero)
    return energy 
end
#=}}}=#


"""
    compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par; thresh=1e-8, max_depth=20)
"""
function compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par; thresh=1e-8, max_depth=20)
    #={{{=#

    energy = 0.0
    npath_nonzero = 0
    npath_zero = 0
    #hi = 2
    #println(ham_par[hi])
    for hi in 1:length(ham_ops)
        f = iterate_dfs!(ref_state, 
                         ham_ops[hi], ham_par[hi], 
                         ansatz_ops, ansatz_par, 
                         thresh=thresh, max_depth=max_depth)

        energy += f[1]
        npath_nonzero += f[2]
        npath_zero += f[3]
    end
    @printf(" E(cADAPT) = %12.8f  #Diagonal Paths %12i  #Nondiagonal Paths %12i\n", energy, npath_nonzero, npath_zero)
    return energy 
end
#=}}}=#


"""
"""
function iterate_dfs!(ref_state, o::P, h::T, ansatz_ops::Vector{Vector{P}}, 
                      ansatz_par::Vector{T}; thresh=1e-12, max_depth=20) where {T,N, P<:Pauli}
#={{{=#
    vcos = cos.(2 .* ansatz_par)
    vsin = sin.(2 .* ansatz_par)
    depth = 0
   
    stack = Stack{Tuple{typeof(o),Float64,Int,Int}}()  
    #stack = Stack{Tuple{typeof(o),Float64,Int,Int}}(undef,1000)  
    #stack = Vector{Tuple{typeof(o),Float64,Int,Int}}()  

    
    push!(stack, (o,h,1,1)) 

    my_energy::Vector{T} = [0.0]
    my_paths::Vector{Int} = [0, 0]

    while length(stack) > 0
        oi, hi, ansatz_layer, depth = pop!(stack)
    
        if ansatz_layer == length(ansatz_ops)+1
            _found_leaf(ref_state, my_energy, oi, hi, my_paths )
        elseif abs(hi) < thresh
            _found_leaf(ref_state, my_energy, oi, hi, my_paths )
        else
#            start = ansatz_layer
#            for i in start:length(ansatz_ops)-1
#                if commute(ansatz_ops[ansatz_layer],oi) 
#                    ansatz_layer += 1
#                else
#                    break 
#                end
#            end
#            
#            g = ansatz_ops[ansatz_layer]
#            if ansatz_layer == length(ansatz_ops)+1
#                _found_leaf(ref_state, my_energy, oi, hi, my_paths )
#            else
#                if depth >= max_depth
#                    _found_leaf(ref_state, my_energy, oi, hi, my_paths )
#                end
#
#                phase, or = commutator(g, oi)
#                hr = real(1im*phase) * hi * vsin[ansatz_layer]
#                #hr = 0.5*real(1im*phase) * hi * vsin[ansatz_layer]
#
#                push!(stack, (or, hr, ansatz_layer+1, depth+1))
#
#                # left branch
#                hl = hi * vcos[ansatz_layer]
#                push!(stack, (oi, hl, ansatz_layer+1, depth))
#            end

            g = ansatz_ops[ansatz_layer]
            if commute(g,oi)
                push!(stack, (oi, hi, ansatz_layer+1, depth))
            else
                if depth >= max_depth
                    _found_leaf(ref_state, my_energy, oi, hi, my_paths )
                end

                phase, or = commutator(g, oi)
                hr = real(1im*phase) * hi * vsin[ansatz_layer]
                #hr = 0.5*real(1im*phase) * hi * vsin[ansatz_layer]

                push!(stack, (or, hr, ansatz_layer+1, depth+1))

                # left branch
                hl = hi * vcos[ansatz_layer]
                push!(stack, (oi, hl, ansatz_layer+1, depth))
            end
        end
    end
    return my_energy[1], my_paths[1], my_paths[2] 
end
#=}}}=#

"""
"""
function iterate_dfs!(ref_state, o::P, h::T, ansatz_ops::Vector{P}, 
                      ansatz_par::Vector{T}; thresh=1e-12, max_depth=20) where {T,N, P<:Pauli}
#={{{=#
    vcos = cos.(2 .* ansatz_par)
    vsin = sin.(2 .* ansatz_par)
    depth = 0
   
    stack = Stack{Tuple{typeof(o),Float64,Int,Int}}()  
    #stack = Stack{Tuple{typeof(o),Float64,Int,Int}}(undef,1000)  
    #stack = Vector{Tuple{typeof(o),Float64,Int,Int}}()  

    
    push!(stack, (o,h,1,1)) 

    my_energy::Vector{T} = [0.0]
    my_paths::Vector{Int} = [0, 0]

    while length(stack) > 0
        oi, hi, ansatz_layer, depth = pop!(stack)
    
        if ansatz_layer == length(ansatz_ops)+1
            _found_leaf(ref_state, my_energy, oi, hi, my_paths )
        elseif abs(hi) < thresh
            _found_leaf(ref_state, my_energy, oi, hi, my_paths )
        else
#            start = ansatz_layer
#            for i in start:length(ansatz_ops)-1
#                if commute(ansatz_ops[ansatz_layer],oi) 
#                    ansatz_layer += 1
#                else
#                    break 
#                end
#            end
#            
#            g = ansatz_ops[ansatz_layer]
#            if ansatz_layer == length(ansatz_ops)+1
#                _found_leaf(ref_state, my_energy, oi, hi, my_paths )
#            else
#                if depth >= max_depth
#                    _found_leaf(ref_state, my_energy, oi, hi, my_paths )
#                end
#
#                phase, or = commutator(g, oi)
#                hr = real(1im*phase) * hi * vsin[ansatz_layer]
#                #hr = 0.5*real(1im*phase) * hi * vsin[ansatz_layer]
#
#                push!(stack, (or, hr, ansatz_layer+1, depth+1))
#
#                # left branch
#                hl = hi * vcos[ansatz_layer]
#                push!(stack, (oi, hl, ansatz_layer+1, depth))
#            end

            g = ansatz_ops[ansatz_layer]
            if commute(g,oi)
                push!(stack, (oi, hi, ansatz_layer+1, depth))
            else
                if depth >= max_depth
                    _found_leaf(ref_state, my_energy, oi, hi, my_paths )
                end

                phase, or = commutator(g, oi)
                hr = real(1im*phase) * hi * vsin[ansatz_layer]
                #hr = 0.5*real(1im*phase) * hi * vsin[ansatz_layer]

                push!(stack, (or, hr, ansatz_layer+1, depth+1))

                # left branch
                hl = hi * vcos[ansatz_layer]
                push!(stack, (oi, hl, ansatz_layer+1, depth))
            end
        end
    end
    return my_energy[1], my_paths[1], my_paths[2] 
end
#=}}}=#


function iterate_dfs2!(ref_state, o::P, h::T, ansatz_ops::Vector{P}, 
                      ansatz_par::Vector{T}; thresh=1e-12, max_depth=20) where {T,N, P<:Pauli}
#={{{=#
    vcos = cos.(2 .* ansatz_par)
    vsin = sin.(2 .* ansatz_par)
    depth = 0
   
    #ori = PauliString(N)
    my_energy::Vector{T} = [0.0]

    stack = Stack{Tuple{typeof(o),Float64,Int}}()  
    push!(stack, (o,h,1)) 

    while length(stack) > 0
        oi, hi, ansatz_layer = pop!(stack)
    
        if ansatz_layer == length(ansatz_ops)+1
            _found_leaf(ref_state, my_energy, oi, hi, paths )
        #elseif abs(hi) < thresh
        #    _found_leaf(ref_state, energy, oi, hi, paths )
        else
            g = ansatz_ops[ansatz_layer]
            if commute(g,oi)
                push!(stack, (oi, hi, ansatz_layer+1))
            else
                #if depth < max_depth
                #if abs(vsin[ansatz_layer]*h) > thresh

                # right branch
                #@btime  commutator!($g, $oi, $or)
                #@code_warntype  commutator(g, oi)
                #error("here")
                
                #phase = commutator!(g, oi, ori)
                #phase2, or2 = commutator(g, oi)
                #phase2 == phase || error(" phase:", phase2, phase)
                #or2 == ori || error(" ori:", g, oi, or2, ori)

                if abs(hi * vsin[ansatz_layer]) > thresh
                    phase, or = commutator(g, oi)
                    hr = real(1im*phase) * hi * vsin[ansatz_layer]
                    #hr = 0.5*real(1im*phase) * hi * vsin[ansatz_layer]

                    push!(stack, (or, hr, ansatz_layer+1))

                    # left branch
                    hl = hi * vcos[ansatz_layer]
                    push!(stack, (oi, hl, ansatz_layer+1))
                else
                    #_found_leaf(ref_state, energy, oi, hi, paths )
                    push!(stack, (oi, hi, ansatz_layer+1))
                end
            end
        end
    end
    return my_energy[1]
end
#=}}}=#


