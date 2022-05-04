using UnitaryPruning
using Plots
using PyCall
using Distributed 
using LinearAlgebra
using Printf

np = pyimport("numpy")
pickle = pyimport("pickle")

T = UInt

dir = "./"
ref_state = [1,1,1,1,1,1,0,0,0,0,0,0,0,0]
N = length(ref_state) 


if nprocs() > 1
    @everywhere ref_state = $ref_state 
    @everywhere include("./read_in.jl")
    @everywhere T = $T
    @everywhere N = $N
else
    include("./read_in.jl")
end

e_list = []
e_list2 = []
thresh_list = []
t_list = []
for i in 2:8
    thresh = 10.0^(-i)
    #t = @elapsed ei = UnitaryPruning.compute_expectation_value_iter_parallel(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh)
    t = @elapsed ei,gi = UnitaryPruning.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, clip=1e-8, thresh=thresh, verbose=1)
    push!(e_list, ei)
    #push!(e_list2, ei2)
    push!(thresh_list, thresh)
    push!(t_list, t)
end
for ti in 1:length(t_list)
    @printf(" %12.8f %12.1e\n", e_list[ti], t_list[ti])
end
e_hf = -17.084018547385366
ref_val = -17.24509797
#ref_val = -2.23298977
#ref_val = -2.23298944

#for i in 2:4
#    thresh = 10.0^(-i)
#    res = UnitaryPruning.optimize_params(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=1e-14, thresh1=thresh);
#end

#plot(thresh_list, [e_list .- ref_val, e_list2 .- ref_val, [e_hf-ref_val for i in 1:length(e_list)]],  
#plot(thresh_list, [e_list .- ref_val, [e_hf-ref_val for i in 1:length(e_list)]],  
plot(thresh_list, [abs.(e_list .- ref_val), ],  
     label = ["classical ADAPT",], 
     lw = 3, 
     marker=true,
     xaxis=:log,
     yaxis=:log)
xlabel!("Threshold")
ylabel!("Error, au")
title!("BeH2 @ 2.39A")
