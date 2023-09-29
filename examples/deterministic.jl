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

function run(; N=6)
   
    ket = KetBitString(N, 0) 
#    o = Pauli(N, Z=[5])
#    o = Pauli(N, Y=[1])
#    o = Pauli(N, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33])
#    o = Pauli(N, Z=[2,5,7,8], Y=[1,3,4])
    o = Pauli(N, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33])

    angles = zeros(Float64,17)
    e = zeros(ComplexF64,17)
    for i in [(i-1)*2 for i in 1:17]
#    for i in 3:3        
        #
        # Uncomment the following to do a serial run
        #
        α = i * π / 32
        generators, parameters = UnitaryPruning.eagle_processor(o, α=α, k=5)
        angles[Int(i/2)+1] = α
        
        # ei , nops = UnitaryPruning.deterministic_pauli_rotations(generators, parameters, o, ket, thres=1e-3)
        ei , nops = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=1e-3)
       
        e[Int(i/2)+1] += ei
        # @printf(" α: %6.4f e: %12.8f+%12.8fi\n", α, real(e[Int(i/2)+1]), imag(e[Int(i/2)+1]))
        @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %i\n", α, real(e[Int(i/2)+1]), imag(e[Int(i/2)+1]), maximum(nops))
        
    end
    
    # plot(angles, real(e))
#    xlabel!("Angles")
#    ylabel!("expectation value")
#    title!{X_{13,29,31}, Y_{9,30}, Z_{8,12,17,28,32}}
#    savefig("plot.pdf")
    return e
end

@time v,e = run(N=127);

function run_eagle(α, o::Pauli{N}; k=5, ket_idx=0) where N
    N == 127 || throw(DimensionMismatch)
    ket = KetBitString(N, ket_idx) 
    generators, parameters = UnitaryPruning.eagle_processor(o, α=α, k=5)
    # expval, n_ops = UnitaryPruning.deterministic_pauli_rotations(generators, parameters, o, ket, thres=1e-3)
    expval, n_ops = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=1e-3)
    return expval, n_ops
end        

# @time expval, n_ops2 = run_eagle(π/4, Pauli(127, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33]))
    