using UnitaryPruning
using Plots

include("../src/read_in.jl")

ref_val = -3.04433127
e_hf = -2.7221671004


ref_state = [1,1,1,1,0,0,0,0]

e_list = []
thresh_list = []
for i in 2:9
    thresh = 10.0^(-i)
    @time ei = UnitaryPruning.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh)
    #@time ei = up.compute_expectation_value_recurse(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=1e-4)
    push!(e_list, ei)
    push!(thresh_list, thresh)
end

plot(thresh_list, [e_list .- ref_val, [e_hf-ref_val for i in 1:length(e_list)]],  
     label = ["classical ADAPT" "HF"], 
     lw = 3, 
     marker=true,
     xaxis=:log)
xlabel!("Threshold")
ylabel!("Error, au")
title!("H4 @ 2A")
