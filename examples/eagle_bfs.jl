using Distributed
@everywhere begin
    using UnitaryPruning
    # using Plots
    using Statistics
    using Printf
    using Random
    using LinearAlgebra
    using SharedArrays
    using PauliOperators

end

function run(; N=6, thresh=1e-3)
   
    ket = KetBitString(N, 0) 
#    o = Pauli(N, Z=[5])
#    o = Pauli(N, Y=[1])
#    o = Pauli(N, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33])
#    o = Pauli(N, Z=[2,5,7,8], Y=[1,3,4])
    o = Pauli(N, Z=[63])
    k = 20 
    
    o = Pauli(N, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33])
    k = 5


    e = Vector{ComplexF64}([]) 
    angles = Vector{Float64}([]) 
    for i in 0:16 
        
        α = i * π / 32
        generators, parameters = UnitaryPruning.eagle_processor(o, α=α, k=k)
        
        # ei , nops = UnitaryPruning.deterministic_pauli_rotations(generators, parameters, o, ket, thres=1e-3)
        ei , nops, l1, l2, l4, entropy = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)
       
        push!(e, ei)
        push!(angles, α)
        
        # @printf(" α: %6.4f e: %12.8f+%12.8fi\n", α, real(e[Int(i/2)+1]), imag(e[Int(i/2)+1]))
        @printf(" α: %6.4f e: %12.8f +%12.8fi l1: %12.8f l2: %12.8f l4: %12.8f shannon: %12.8f nops: %i\n", α, real(ei), imag(ei), l1, l2, 1-l4, entropy, maximum(nops))
        
    end
    
    plot(angles, real(e), marker=:circle)
    xlabel!("Angles")
    ylabel!("expectation value")
#    title!{X_{13,29,31}, Y_{9,30}, Z_{8,12,17,28,32}}
    savefig("plot_eagle_bfs.pdf")
    return e
end

@time v,e = run(N=127, thresh=1e-2);

function run_eagle(α, o::Pauli{N}; k=5, ket_idx=0) where N
    N == 127 || throw(DimensionMismatch)
    ket = KetBitString(N, ket_idx) 
    generators, parameters = UnitaryPruning.eagle_processor(o, α=α, k=5)
    # expval, n_ops = UnitaryPruning.deterministic_pauli_rotations(generators, parameters, o, ket, thres=1e-3)
    expval, n_ops = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=1e-3)
    return expval, n_ops
end        

# @time expval, n_ops2 = run_eagle(π/4, Pauli(127, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33]))
    