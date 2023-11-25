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

function run(; N=6, k=10, thresh=1e-3)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    angles = [] 
    e = [] 
    
    # for i in [(i-1)*2 for i in 1:9]
    for i in 0:16
    # for i in 8:8
        α = i * π / 32 
        generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=k)
        
        ei , nops, l1, l2, l4, entropy = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)
       
        push!(e, ei)
        push!(angles, α)
        
        # @printf(" α: %6.4f e: %12.8f+%12.8fi\n", α, real(e[Int(i/2)+1]), imag(e[Int(i/2)+1]))
        @printf(" α: %6.4f e: %12.8f +%12.8fi l1: %12.8f l2: %12.8f l4: %12.8f shannon: %12.8f nops: %i\n", α, real(ei), imag(ei), l1, l2, 1-l4, entropy, maximum(nops))
        
    end
    
    plot(angles, real(e), marker = :circle)
    xlabel!("Angle")
#    ylabel!("expectation value")
#    title!{X_{13,29,31}, Y_{9,30}, Z_{8,12,17,28,32}}
    savefig("plot_1d_n6_bfs.pdf")
    return e
end

@time v,e = run(k=5, N=6, thresh=1e-4);
