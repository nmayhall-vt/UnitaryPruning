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

function run(;thresh=1e-3)
  
    N = 127

    ket = KetBitString(N, 0) 
    o = Pauli(N, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33])

    e = Vector{ComplexF64}([]) 
    angles = Vector{Float64}([]) 
    # for i in 0:16 
    for i in 8:8 
        α = i * π / 32
        generators, parameters = UnitaryPruning.eagle_processor(o, α=α, k=5)
        ei, l2a, l2b = UnitaryPruning.dfs_evolution(generators, parameters, 1.0, o, ket, thresh=thresh)
        
        push!(e, ei)
        push!(angles, α)
        @printf(" α: %6.4f e: %12.8f +%12.8fi l2a: %12.8f l2b: %12.8f l2: %12.8f\n", α, real(ei), imag(ei), l2a, l2b, sqrt(l2a+l2b))
    end
    
    # plot(angles, real(e), marker=:circle)
    # xlabel!("Angles")
    # ylabel!("expectation value")
    # # title!{X_{13,29,31}, Y_{9,30}, Z_{8,12,17,28,32}}
    # savefig("plot_eagle_dfs.pdf")
    return e
end

@time e = run(thresh=1e-4)
