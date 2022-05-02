using UnitaryPruning
using Plots
using PyCall
using Distributed 
using LinearAlgebra
using Printf

np = pyimport("numpy")
pickle = pyimport("pickle")

@everywhere T = UInt

#dir = "./examples/beh2_2r/"
dir = "./"
@everywhere ref_state = [1,1,1,1,1,1,0,0,0,0,0,0,0,0]
@everywhere N = length(ref_state) 
include("./read_in.jl")

@everywhere ham_ops = $ham_ops
@everywhere ham_par = $ham_par
@everywhere ansatz_ops = $ansatz_ops
@everywhere ansatz_par = $ansatz_par


#function mypickle(filename, obj)
#    out = open(filename,"w")
#    pickle.dump(obj, out)
#    close(out)
# end
#
#function myunpickle(filename)
#    r = nothing
#    @pywith pybuiltin("open")(filename,"rb") as f begin
#        r = pickle.load(f)
#    end
#    return r
#end
#
#grouped = myunpickle(dir*"ansatz_ops_grouped.pkl")
#for i in grouped
#    tmp::Vector{Tuple{PauliString{N},Float64}} = []
#    for j in i
#        println(j[2], " ", abs(j[1]))
#        push!(tmp, (PauliString(j[2]), real(-1.0im*j[1])))
#    end
#    push!(ansatz_ops_grouped, tmp)
#    println()
#end
#
#for a in ansatz_ops_grouped
#    for i in a
#        for j in a
#            commute(i[1], j[1]) || error("here")
#        end
#    end
#end


e_list = []
e_list2 = []
thresh_list = []
t_list = []
for i in 1:6
    thresh = 10.0^(-i)
    #t = @elapsed ei = UnitaryPruning.compute_expectation_value_iter_parallel(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh)
    t = @elapsed ei,gi = UnitaryPruning.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=1e-8, thresh1=thresh)
    #t = @elapsed ei = UnitaryPruning.compute_expectation_value_recurse(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh)
    println(ei)
    push!(e_list, ei)
    #push!(e_list2, ei2)
    push!(thresh_list, thresh)
    push!(t_list, t)
end
for ti in 1:length(t_list)
    @printf(" %12.8f %12.1e\n", e_list[ti], t_list[ti])
end
e_hf = -17.084018547385366
ref_val = -17.21643699


#for i in 2:4
#    thresh = 10.0^(-i)
#    res = UnitaryPruning.optimize_params(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=1e-14, thresh1=thresh);
#end

#plot(thresh_list, [e_list .- ref_val, e_list2 .- ref_val, [e_hf-ref_val for i in 1:length(e_list)]],  
plot(thresh_list, [e_list .- ref_val, [e_hf-ref_val for i in 1:length(e_list)]],  
     label = ["classical ADAPT" "HF"], 
     lw = 3, 
     marker=true,
     xaxis=:log)
xlabel!("Threshold")
ylabel!("Error, au")
title!("BeH2 @ 2.39A")
