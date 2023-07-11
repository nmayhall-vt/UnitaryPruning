@everywhere begin
    using Distributed
    using UnitaryPruning
    using Plots
    using Statistics
    using Printf
    using Random
    using LinearAlgebra
    using SharedArrays
end

@everywhere function run(;α=.01, k=10, nsamples=10000, seed=1)
    
    N = 6  # number of qubits

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

    # o = PauliBoolVec(N,X=[1],Y=[2],Z=[3])

    #Mz
    # o = PauliBoolVec(N, X=[13,29,31], Y=[9,30], Z=[8,12,17,28,32])
    # o = PauliBoolVec(N, Z=[1,2,3,4])
    o = PauliBoolVec(N, Y=[1], Z=[2,3,4])
    # o = PauliBoolVec(N, Z=[1,30])
    # o = PauliBoolVec(N, Z=[1])


    ket = zeros(Bool, N)
    Random.seed!(seed)
    rolling_avg = zeros(ComplexF64, nsamples)
    
    samples = UnitaryPruning.stochastic_pauli_dynamics(generators, parameters, o, ket, nsamples=nsamples)

    rolling_avg[1] = samples[1] 
    
    for i in 2:nsamples
        rolling_avg[i] = (rolling_avg[i-1]*(i-1) + samples[i])/i
    end
    return rolling_avg, samples 
end



nruns = 100
nsamples = 10000

final_vals_stoc = []
final_vals_errs = []

# for i in 7:8
for i in 0:16
    avg_traj = zeros(nsamples)
   
    @everywhere nsamples = $nsamples
    @everywhere trajectory = zeros(nsamples)
    # trajectories = [zeros(nsamples) for i in nruns]


    # trajectories = SharedArray{Vector{Float64}}([zeros(nsamples) for i in nruns])
    # trajectories = SharedMatrix{Float64}(nruns,nsamples)
    # avg_traj = SharedVector{Float64}(nsamples)

    avg_traj = @sync @distributed (.+) for runi in 1:nruns
    # for runi in 1:nruns
        ei, _ = run(α=i*π/32, k=5, nsamples=nsamples, seed=runi)
        ei .= real.(ei)
        traj = ei
        # trajectory .= ei
        # trajectories[runi,:] .= ei
        # avg_traj .+= ei
    end
    avg_traj ./= nruns
    var_traj = sum([(i .- avg_traj).*(i .- avg_traj) for i in trajectories]) ./nruns
    
    std_traj = sqrt.(var_traj)
    z = 1.96 # 95% confidence
    std_traj = z .* std_traj ./ sqrt(nruns)

    @printf(" α: %4i avg: %12.8f ± %-12.6f var: %12.8f\n", i, avg_traj[end], std_traj[end], var_traj[end])

    # plot(trajectories, ylim=[-1,1], legend=false, color=:grey, alpha=.1)
    # plot!(avg_traj, color=:black, ribbon=std_traj)
    # plot!(avg_traj .+ std_traj, color=:red)
    # plot!(avg_traj .- std_traj, color=:red)
    # savefig(@sprintf "%05i.pdf" i)
    push!(final_vals_stoc, avg_traj[end])
    push!(final_vals_errs, std_traj[end])
end

# plot(final_vals_stoc, ribbon=final_vals_errs, marker=:circle); 
# plot!(final_vals); 
# savefig("plot.pdf")