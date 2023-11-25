
@everywhere function compute_run(generators, parameters, o, ket, sample_size; seed=1, verbose=1)
    Random.seed!(seed)
    
    sample = UnitaryPruning.stochastic_pauli_rotations(generators, parameters, o, ket, sample_size=sample_size)

    expval = 0.0
    l2 = 0.0
    entropy = 0.0
    for (oi,count) in sample 
        prob_i = count / sample_size
        
        entropy -= log2(prob_i)*prob_i
        
        expval += sqrt(prob_i)*expectation_value(oi, ket)
        
        l2 += prob_i
    end
   
    verbose < 1 || @printf(" ev: %12.8f entropy: %12.8f #ops: %5i l2: %12.8f\n", expval, entropy, length(sample), l2)

    return expval
    # return expval, entropy
    # # @printf(" ev: %12.8f l2: %12.8f normalized <o>: %12.8f size of pop: %5i\n", expval, l2, expval/l2, length(population))
    
    # rolling_avg[1] = samples[1]
    
    # for i in 2:sample_size
    #     rolling_avg[i] = (rolling_avg[i-1]*(i-1) + samples[i])/i
    # end
    # return rolling_avg, rolling_avg.*rolling_avg 
end


function run(; nruns=100, sample_size=1000, k=5, ket_idx=0)
    N = 127 
    o = Pauli(N, Z=[1])
    o = Pauli(N, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33])
    
    ket = KetBitString(N, ket_idx) 
   

    final_vals_stoc = []
    final_vals_errs = []

    # for i in [(i-1)*2 for i in 1:9]
    # for i in 0:16
    for i in 8:8 
    
        α = i * π / 32
        generators, parameters = UnitaryPruning.eagle_processor(o, α=α, k=k)

        avg_vec = zeros(nruns)
        var_vec = zeros(nruns)
        std_vec = zeros(nruns)
        samples = zeros(nruns)
        for runi in 1:nruns
            samples[runi] = compute_run(generators, parameters, o, ket, sample_size, seed=runi)
            avg_vec[runi] = mean(samples[1:runi]) 
            var_vec[runi] = var(samples[1:runi])
            std_vec[runi] = std(samples[1:runi])
        end

        
        std_vec = sqrt.(abs.(var_vec))
        z = 1.96 # 95% confidence
        std_vec = z .* std_vec ./ sqrt(nruns)

        @printf(" α: %12.8f avg: %12.8f ± %-12.6f var: %12.8f\n", α, mean(samples), std(samples)*z, var(samples))

        #
        # Plot stuff
        #
        # plot(vecectories, ylim=[-1,1], legend=false, color=:grey, alpha=.1)
        # plot([mean(avg_vec[1:i]) for i in 1:length(avg_vec)], color=:black, ribbon=std_vec, ylim=[-1, 1], legend=false, dpi=300)
        plot([mean(avg_vec[1:i]) for i in 1:length(avg_vec)], color=:black, ribbon=std_vec, legend=false, dpi=300)
        filename = @sprintf "%05i.png" i
        savefig(filename)
        push!(final_vals_stoc, avg_vec[end])
        push!(final_vals_errs, std_vec[end])
    end

    plot(final_vals_stoc, ribbon=final_vals_errs, marker=:circle)
    # plot!(final_vals)
    savefig("plot.pdf")
    return final_vals_stoc, final_vals_errs
end


@time v,e = run(nruns=10, sample_size=100000000)

