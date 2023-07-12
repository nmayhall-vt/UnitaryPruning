using Distributed
@everywhere begin
    using UnitaryPruning
    using Plots
    using Statistics
    using Printf
    using Random
    using LinearAlgebra
    using SharedArrays
end

@everywhere function get_unitary_sequence_1D(;α=.01, k=10, N=6)
    
    generators = Vector{PauliBoolVec{N}}([])
    parameters = Vector{Float64}([])

    # Loop over trotter steps
    for ki in 1:k
        
        ## ZZ layer
        # e^{i π/2 P2} e^{i π P1 /2}|ψ>
        for i in 1:N-1
            pi = PauliBoolVec(N, Z=[i, i + 1])
            push!(generators, pi)
            push!(parameters, π)
        end
            
        #pbc 
        pi = PauliBoolVec(N, Z=[N, 1])
        push!(generators, pi)
        push!(parameters, π )

        ## X layer
        # e^{i αn Pn / 2}
        for i in 1:N
            pi = PauliBoolVec(N, X=[i])
            push!(generators, pi)
            push!(parameters, α)
        end
    end

    return generators, parameters
end

@everywhere function generate_samples(generators, parameters, o, ket, nsamples; seed=1)
    Random.seed!(seed)
    rolling_avg = zeros(nsamples)
    
    samples = UnitaryPruning.randomized_pauli_dynamics(generators, parameters, o, ket, nsamples=nsamples)
    samples .= real.(samples)

    rolling_avg[1] = samples[1]
    
    for i in 2:nsamples
        rolling_avg[i] = (rolling_avg[i-1]*(i-1) + samples[i])/i
    end
    return rolling_avg, rolling_avg.*rolling_avg 
end


function run(;N=100, nruns=100, nsamples=1000)

    final_vals_stoc = []
    final_vals_errs = []

    for i in 8:8

        # Operator
        # o = PauliBoolVec(N, X=[13,29,31], Y=[9,30], Z=[8,12,17,28,32])
        # o = PauliBoolVec(N, Z=[1,2,3,4])
        o = PauliBoolVec(N, Y=[1], Z=[2,3,4])
        # o = PauliBoolVec(N, Z=[1,30])
        o = PauliBoolVec(N, Z=[1,2])
        o = PauliBoolVec(N, Z=[1,2,3,4,5,6])

        # State
        ket = zeros(Bool, N)

        avg_traj = zeros(nsamples)
        var_traj = zeros(nsamples)
        std_traj = zeros(nsamples)

        @everywhere generators, parameters = get_unitary_sequence_1D(α=$i * π / 32, k=5, N=$N)

        avg_traj, var_traj = @sync @distributed (.+) for runi in 1:nruns
            generate_samples(generators, parameters, o, ket, nsamples, seed=runi)
        end

        # a = pmap(runi -> generate_samples(generators, parameters, o, ket, nsamples, seed=runi), 1:nruns)

        # for ai in fetch(a)
        #     # println(size(ai[1]), size(ai[2]))
        #     avg_traj .+= ai[1]
        #     var_traj .+= ai[1].*ai[1]
        # end

        var_traj .= var_traj .- (avg_traj .* avg_traj) ./ nruns
        avg_traj ./= nruns
        var_traj ./= nruns
        # var_traj = sum([(i .- avg_traj).*(i .- avg_traj) for i in trajectories]) ./nruns

        std_traj = sqrt.(var_traj)
        z = 1.96 # 95% confidence
        std_traj = z .* std_traj ./ sqrt(nruns)

        @printf(" α: %4i avg: %12.8f ± %-12.6f var: %12.8f\n", i, avg_traj[end], std_traj[end], var_traj[end])

        # plot(trajectories, ylim=[-1,1], legend=false, color=:grey, alpha=.1)
        plot(avg_traj, color=:black, ribbon=std_traj, ylim=[-1, 1], legend=false, dpi=300)
        # plot!(avg_traj .+ std_traj, color=:red)
        # plot!(avg_traj .- std_traj, color=:red)
        filename = @sprintf "%05i.png" i
        savefig(filename)
        push!(final_vals_stoc, avg_traj[end])
        push!(final_vals_errs, std_traj[end])
    end

    plot(final_vals_stoc, ribbon=final_vals_errs, marker=:circle)
    plot!(final_vals)
    savefig("plot.pdf")
end




run(N=6, nruns=1000, nsamples=10000000)