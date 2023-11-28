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



function compute_run(generators, parameters, o, ket, sample_size; seed=1, verbose=1)
    Random.seed!(seed)
    
    sample, H = UnitaryPruning.stochastic_pauli_rotations(generators, parameters, o, ket, sample_size=sample_size)

    expval = 0.0
    l2 = 0.0
    entropy = 0.0
    for (oi,count) in sample 
        prob_i = count / sample_size
        
        entropy -= log2(prob_i)*prob_i
        
        expval += sqrt(prob_i) * expectation_value(oi, ket)
        
        l2 += prob_i
    end
 
    abs(l2-1.0) < 1e-14 || @warn "normalization: ", l2
    verbose < 1 || @printf(" ev: %12.8f entropy: %12.8f #ops: %5i \n", expval, entropy, length(sample))

    return expval
end


function run(; nruns=100, sample_size=1000, N=6, k=5)
   
    #
    # Uncomment to use bitstrings
    #
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    o_mat = Matrix(o)

    e_ref = []
    final_vals_stoc = []
    final_vals_errs = []

    # for i in [(i-1)*2 for i in 1:9]
    # for i in 0:16
    # for i in 15:15
    for i in 8:8
        α = i*π/32
        generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=k)
       
        g2 = similar(generators)
        p2 = similar(parameters) 
        n_repeat = 1 
        for g in generators
            for i in 1:n_repeat
                push!(g2, g)
            end
        end
        for p in parameters 
            for i in 1:n_repeat
                push!(p2, p/n_repeat)
            end
        end
        generators = g2
        parameters = p2


        U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
        ei = diag(U'*o_mat*U)[1]

        push!(e_ref, real(ei))
        
        avg_vec = zeros(nruns)
        var_vec = zeros(nruns)
        std_vec = zeros(nruns)
        samples = zeros(nruns)
        conf_interval = zeros(nruns)
        
        z = 1.96 # 95% confidence
        
        for runi in 1:nruns
            samples[runi] = compute_run(generators, parameters, o, ket, sample_size, seed=runi)
            avg_vec[runi] = mean(samples[1:runi]) 
            var_vec[runi] = var(samples[1:runi])
            std_vec[runi] = std(samples[1:runi])
            conf_interval[runi] = z .* std_vec[runi] ./ sqrt(runi)
        end

        
        @printf(" α: %12.8f avg: %12.8f ± %-12.6f var: %12.8f Eref: %12.8f\n", α, mean(samples), std(samples)*z/sqrt(nruns), var(samples), ei)


        #
        # Plot stuff
        #
        # plot(trajectories, ylim=[-1,1], legend=false, color=:grey, alpha=.1)
        # plot(avg_traj, color=:black, ribbon=std_traj, ylim=[-1, 1], legend=false, dpi=300)
        # filename = @sprintf "%05i.png" i
        # savefig(filename)
        push!(final_vals_stoc, avg_vec[end])
        push!(final_vals_errs, std_vec[end])
    end

    plot(final_vals_stoc, ribbon=final_vals_errs, marker=:circle)
    xlabel!("Angle")
    savefig("plot_1d_n6_stochastic.pdf")
    return final_vals_stoc, final_vals_errs
end


@time v,e = run(nruns=2, sample_size=100000, k=5)
