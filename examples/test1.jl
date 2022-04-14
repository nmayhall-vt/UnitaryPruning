using UnitaryPruning
using Plots
using PyCall

np = pyimport("numpy")

T = UInt
N = 12

println(dir)

ham_par    = Vector{Float64}()
ansatz_par = Vector{Float64}()

ham_ops    = Vector{PauliString{N}}()
ansatz_ops = Vector{PauliString{N}}()

ref_state = [1,1,1,1,1,1,0,0,0,0,0,0]


for i in np.load(dir*"/ansatz_par.npy")
    # Ansatz parameters should be real 
    # openfermion tends to add phase to coefficients
    push!(ansatz_par, real(-1im * i))
end
for i in np.load(dir*"/ham_par.npy")
    push!(ham_par, i)
end

do_mask = true 

if do_mask
    ansatz_ops = Vector{PauliMask{T,N}}()
    ham_ops    = Vector{PauliMask{T,N}}()
    ref_state = parse(T, join(ref_state); base=2)
    for i in np.load(dir*"/ham_ops.npy")
        push!(ham_ops, PauliMask(PauliString(i)))
    end
    for i in np.load(dir*"/ansatz_ops.npy")
        push!(ansatz_ops, PauliMask(PauliString(i)))
    end
else
    for i in np.load(dir*"/ham_ops.npy")
        push!(ham_ops, PauliString(i))
    end
    for i in np.load(dir*"/ansatz_ops.npy")
        push!(ansatz_ops, PauliString(i))
    end
end

ref_val = -3.04433127
e_hf = -2.7221671004

ref_val = -7.83360160
e_hf = -7.7393739490



e_list = []
e_list2 = []
thresh_list = []
t_list = []
for i in 1:9
    thresh = 10.0^(-i)
    t = @elapsed ei = UnitaryPruning.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=thresh)
    push!(e_list, ei)
    #push!(e_list2, ei2)
    push!(thresh_list, thresh)
    push!(t_list, t)
end
for t in t_list
    println(t)
end

#plot(thresh_list, [e_list .- ref_val, e_list2 .- ref_val, [e_hf-ref_val for i in 1:length(e_list)]],  
plot(thresh_list, [e_list .- ref_val, [e_hf-ref_val for i in 1:length(e_list)]],  
     label = ["classical ADAPT" "HF"], 
     lw = 3, 
     marker=true,
     xaxis=:log)
xlabel!("Threshold")
ylabel!("Error, au")
title!("H6 @ 1A")
