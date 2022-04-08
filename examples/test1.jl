using UnitaryPruning
using Plots

include("../src/read_in.jl")
e_list = []
thresh_list = []
for i in 1:10
    thresh = 10.0^(-i)
    @time ei = up.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh)
    #@time ei = up.compute_expectation_value_recurse(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=1e-4)
    push!(e_list, ei)
    push!(thresh_list, thresh)
end

plot(thresh_list, [e_list .- ref_val, [-4.3916471843+4.45948534 for i in 1:length(e_list)]],  
     label = ["classical ADAPT" "HF"], 
     lw = 3, 
     marker=true,
     xaxis=:log)
xlabel!("Threshold")
ylabel!("Error, au")
title!("H4 @ 1A")
