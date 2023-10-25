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

    angles = zeros(Float64,17)
    e = zeros(ComplexF64,17)
    for i in [(i-1)*2 for i in 1:9]
    # for i in 3:3        
        #
        # Uncomment the following to do a serial run
        #
        α = i * π / 32
        generators, parameters = UnitaryPruning.eagle_processor(o, α=α, k=5)
        angles[i+1] = α
        
        ei, _ = UnitaryPruning.dfs_evolution(generators, parameters, 1.0, o, ket, thresh=thresh)
        e[i+1] = ei
        @printf(" α: %6.4f e: %12.8f+%12.8fi\n", α, real(ei), imag(ei))
        
    end
    
    # plot(angles, real(e))
    # xlabel!("Angles")
    # ylabel!("expectation value")
    # title!{X_{13,29,31}, Y_{9,30}, Z_{8,12,17,28,32}}
    # savefig("plot_eagle_dfs.pdf")
    return e
end

@time v,e = run(thresh=1e-3)
