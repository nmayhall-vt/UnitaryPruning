using UnitaryPruning
using Plots
using PyCall
using Distributed 
using LinearAlgebra

np = pyimport("numpy")
pickle = pyimport("pickle")

@everywhere T = UInt

#dir = "./examples/beh2_2r/"
dir = "./"
@everywhere ref_state = [1,1,1,1,0,0,0,0]
@everywhere N = length(ref_state) 
include("./read_in.jl")

@everywhere ham_ops = $ham_ops
@everywhere ham_par = $ham_par
@everywhere ansatz_ops = $ansatz_ops
@everywhere ansatz_par = $ansatz_par




g_list = []
e_list = []
e_list2 = []
thresh_list = []
t_list = []
for i in 2:9
    thresh = 10.0^(-i)
    #t = @elapsed ei = UnitaryPruning.compute_expectation_value_iter_parallel(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh)
    t = @elapsed ei,gi = UnitaryPruning.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh, verbose=1)
    #t = @elapsed ei = UnitaryPruning.compute_expectation_value_recurse(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh)
    push!(e_list, ei)
    push!(g_list, gi)
    #push!(e_list2, ei2)
    push!(thresh_list, thresh)
    push!(t_list, t)
end
for t in t_list
    println(t)
end
e_hf = -4.5599680231
ref_val = -4.60060768
do_log = false 
if do_log
    plot(thresh_list, [abs.(e_list .- ref_val), [e_hf-ref_val for i in 1:length(e_list)]],  
         label = ["classical ADAPT" "HF"], 
         lw = 3, 
         marker=true,
         xaxis=:log, yaxis=:log)
else
    plot(thresh_list, [e_list .- ref_val, [e_hf-ref_val for i in 1:length(e_list)]],  
         label = ["classical ADAPT" "HF"], 
         lw = 3, 
         marker=true,
         xaxis=:log)
end
xlabel!("Threshold")
ylabel!("Error, au")
title!("H4 @ 1A")



if false
    e_list = []
    thresh_list = []
    println(" Now we will optimize the parameters for each truncation level")
    for i in 2:8
        thresh = 10.0^(-i)
        res = UnitaryPruning.optimize_params(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, clip=1e-14, thresh=thresh);
        #display(res.minimizer)
        push!(e_list, res.minimum)
        push!(thresh_list, thresh)
    end


    if do_log
        plot!(thresh_list, [abs.(e_list .- ref_val)],  
              label = ["classical ADAPT opt"], 
              lw = 3, 
              marker=true,
              xaxis=:log, yaxis=:log)
    else
        plot!(thresh_list, [e_list .- ref_val],  
              label = ["classical ADAPT opt"], 
              lw = 3, 
              marker=true,
              xaxis=:log)
    end
end
