using Distributed 
using LinearAlgebra
using Optim
using SpecialFunctions
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
function compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par; thresh=1e-8, thresh1=1e-6, max_depth=20, verbose=0)
    #={{{=#

    energy = 0.0
    npath_nonzero = 0
    npath_zero = 0
    gradient = deepcopy(ansatz_par)
    fill!(gradient, zero(0))
    #hi = 2
    #println(ham_par[hi])
    for hi in 1:length(ham_ops)
        f = iterate_dfs!(ref_state, 
                         ham_ops[hi], ham_par[hi], 
                         ansatz_ops, ansatz_par, 
                         thresh=thresh, 
                         thresh1=thresh1,
                         max_depth=max_depth)

        energy += f[1]
        gradient .+= f[2]
        npath_nonzero += f[3]
        npath_zero += f[4]
    end
    if verbose>0
        @printf(" E(cADAPT) = %12.8f  Grad = %8.1e #Diagonal Paths %12i  #Nondiagonal Paths %12i\n", energy, norm(gradient), npath_nonzero, npath_zero)
    end
    return energy, gradient 
end
#=}}}=#



"""
"""
function iterate_dfs!(ref_state, o::P, h::T, ansatz_ops::Vector{P}, 
                      ansatz_par::Vector{T}; 
                      thresh=1e-15, 
                      thresh1=1e-12,
                      max_depth=20) where {T,N, P<:Pauli}
#={{{=#
    vcos = cos.(2 .* ansatz_par)
    vsin = sin.(2 .* ansatz_par)
    vcot = cot.(2 .* ansatz_par)
    vtan = tan.(2 .* ansatz_par)
    depth = 0

    path = zeros(Int8, length(ansatz_ops)+1)
   
    stack = Stack{Tuple{typeof(o),Float64,Int,Int, Int8}}()  
    #stack = Stack{Tuple{typeof(o),Float64,Int,Int}}(undef,1000)  
    #stack = Vector{Tuple{typeof(o),Float64,Int,Int}}()  

    
    push!(stack, (o,h,1,1,0)) 

    my_energy::Vector{T} = [0.0]
    my_paths::Vector{Int} = [0, 0]
    
    my_gradient = deepcopy(ansatz_par)
    fill!(my_gradient, zero(0))

    while length(stack) > 0
        oi, hi, ansatz_layer, depth, branch_direction = pop!(stack)
        path[ansatz_layer] = branch_direction

        if ansatz_layer == length(ansatz_ops)+1
            _found_leaf(ref_state, my_energy, my_gradient, oi, hi, path, my_paths, vtan, vcot, thresh1)

        #elseif abs(hi) < thresh
        elseif abs(hi * erf(hi*hi/thresh1/thresh1)) < thresh
            _found_leaf(ref_state, my_energy, my_gradient, oi, hi, path, my_paths, vtan, vcot, thresh1)
        else

            g = ansatz_ops[ansatz_layer]
            if commute(g,oi)
                push!(stack, (oi, hi, ansatz_layer+1, depth, 0))
            else
                if depth >= max_depth
                    _found_leaf(ref_state, my_energy, my_gradient, oi, hi, path, my_paths, vtan, vcot, thresh1)
                end

                phase, or = commutator(g, oi)
                hr = real(1im*phase) * hi * vsin[ansatz_layer]
                #hr = 0.5*real(1im*phase) * hi * vsin[ansatz_layer]

                push!(stack, (or, hr, ansatz_layer+1, depth+1, -1))

                # left branch
                hl = hi * vcos[ansatz_layer]
                push!(stack, (oi, hl, ansatz_layer+1, depth, 1))
            end
        end
    end
    return my_energy[1], my_gradient, my_paths[1], my_paths[2] 
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


function _found_leaf(ref_state, energy::Vector{T}, gradient::Vector{T}, o, h, path::Vector{Int8}, npaths::Vector{Int}, vtan, vcot, thresh1) where {T}
#={{{=#
    if is_diagonal(o)
        sign = expectation_value_sign(o, ref_state) 

        #@printf(" Found energy contribution %12.8f at ansatz layer %5i and depth %5i\n", sign*h, ansatz_layer, depth)
        ei = sign*h
        
        erfe = erf(ei*ei/thresh1/thresh1)
        energy[1] += ei * erfe
        npaths[1] += 1
                
        #dfde = erfe 
        dfde = erfe + 4*ei*ei*exp(-(ei/thresh1)^4)/thresh1/thresh1/sqrt(pi)

        #println(path)
        # compute gradient contribution
        for p in 1:length(vtan)
            if path[p+1] == Int8(1)
                dedx = -2 * ei * vtan[p]
                gradient[p] += dfde * dedx  
                #gradient[p] += -2 * ei * vtan[p]
            elseif path[p+1] == Int8(-1)
                dedx =  2 * ei * vcot[p]
                gradient[p] += dfde * dedx  
                #gradient[p] +=  2 * ei * vcot[p]
            end
        end

    else
        npaths[2] += 1
    end
end
#=}}}=#


function optimize_params(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par; thresh=1e-8, max_depth=20, thresh1=1e-6)

    #ecurr = 0.0
    #gcurr = zeros(length(ansatz_par)) 
    iter = 0
    ecurr, gcurr = UnitaryPruning.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh, thresh1=thresh1)
    
    function func(p::Vector{Float64})
        ecurr, gcurr = UnitaryPruning.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, p, thresh=thresh, thresh1=thresh1)
        return ecurr 
    end
    function grad(g::Vector{Float64}, p::Vector{Float64})
        et, gt = UnitaryPruning.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, p, thresh=thresh, thresh1=thresh1)
        g .= gt
        return gcurr
    end
    function callback(p)
        iter += 1
        @printf(" %5i E = %12.8f G = %12.8f \n", iter, ecurr, norm(gcurr))
        return false 
    end
   
#    alpha = .01
#    p = deepcopy(ansatz_par)
#    for i in 1:100
#        p .-= alpha .* grad(p)
#        func(p)
#    end

    method = "bfgs"    
    gconv = 1e-7
    max_iter = 50 

    options = Optim.Options(
        callback = callback, 
        g_tol=gconv,
        iterations=max_iter,
        store_trace=true, 
        #show_trace=true
    )

    p = deepcopy(ansatz_par)
        
    optmethod = BFGS()
    optmethod = Newton()
    optmethod = ConjugateGradient()
    optmethod = LBFGS()
   
    #res = Optim.optimize(func, p, optmethod, options)
    res = Optim.optimize(func, grad, p, optmethod, options)
    summary(res)
    e = Optim.minimum(res)
    display(res)
    return res
end
