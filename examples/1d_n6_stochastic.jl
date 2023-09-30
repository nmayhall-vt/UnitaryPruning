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



@everywhere function compute_run(generators, parameters, o, ket, nsamples; seed=1)
    Random.seed!(seed)
    rolling_avg = zeros(nsamples)
    
    samples = UnitaryPruning.stochastic_pauli_rotations(generators, parameters, o, ket, nsamples=nsamples)
    samples .= real.(samples)

    rolling_avg[1] = samples[1]
    
    for i in 2:nsamples
        rolling_avg[i] = (rolling_avg[i-1]*(i-1) + samples[i])/i
    end
    return rolling_avg, rolling_avg.*rolling_avg 
end


function run(; nruns=100, nsamples=1000, N=6, k=5)
   
    #
    # Uncomment to use bitstrings
    #
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    final_vals_stoc = []
    final_vals_errs = []

    # for i in [(i-1)*2 for i in 1:9]
    for i in 0:16

        avg_traj = zeros(nsamples)
        var_traj = zeros(nsamples)
        std_traj = zeros(nsamples)

        generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=i * π / 32, k=k)
        for runi in 1:nruns
             a, b = compute_run(generators, parameters, o, ket, nsamples, seed=runi)
             avg_traj .+= a
             var_traj .+= b
        end


        var_traj .= var_traj .- (avg_traj .* avg_traj) ./ nruns
        avg_traj ./= nruns
        var_traj ./= nruns

        std_traj = sqrt.(var_traj)
        z = 1.96 # 95% confidence
        std_traj = z .* std_traj ./ sqrt(nruns)

        @printf(" α: %4i avg: %12.8f ± %-12.6f var: %12.8f\n", i, avg_traj[end], std_traj[end], var_traj[end])

        #
        # Plot stuff
        #
        # plot(trajectories, ylim=[-1,1], legend=false, color=:grey, alpha=.1)
        # plot(avg_traj, color=:black, ribbon=std_traj, ylim=[-1, 1], legend=false, dpi=300)
        # filename = @sprintf "%05i.png" i
        # savefig(filename)
        push!(final_vals_stoc, avg_traj[end])
        push!(final_vals_errs, std_traj[end])
    end

    plot(final_vals_stoc, ribbon=final_vals_errs, marker=:circle)
    xlabel!("Angle")
    savefig("plot_1d_n6_stochastic.pdf")
    return final_vals_stoc, final_vals_errs
end


@time v,e = run(nruns=100, nsamples=1000, k=80)
