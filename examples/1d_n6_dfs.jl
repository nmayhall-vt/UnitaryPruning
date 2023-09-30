using Distributed
@everywhere begin
    using UnitaryPruning
    using Plots
    using Statistics
    using Printf
    using Random
    using LinearAlgebra
    using SharedArrays
    using PauliOperators

end

function run(;N=6, k=5, thresh=1e-3)
  
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])
    
    angles = [] 
    e = [] 

    for i in 0:16
        α = i * π / 32
        generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=k)
        
        ei, _ = UnitaryPruning.dfs_evolution(generators, parameters, 1.0, o, ket, thresh=thresh)
        
        push!(e, ei)
        push!(angles, α)
        
        @printf(" α: %6.4f e: %12.8f+%12.8fi\n", α, real(ei), imag(ei))
        
    end
    
    plot(angles, real(e), marker = :circle)
    xlabel!("Angle")
    savefig("plot_1d_n6_dfs.pdf")
    return e
end

@time v,e = run(N=6, k=10, thresh=1e-4)
