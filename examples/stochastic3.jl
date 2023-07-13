using UnitaryPruning
using Plots
using Statistics
using Printf
using Random
using LinearAlgebra


function eagle_processor(N=127, k=10, α=.01, sequences=[[1,13], [19,32],[38,51],[57,70],[76,89],[95,108],[114,126]],
                       bridges=[[1,15,19],[5,16,23],[9,17,27],[13,18,31],[33,37,52],[29,36,48],[25,35,44],[21,34,40],
                                [38,53,57],[42,54,61],[46,55,65],[50,56,69],[71,75,90],[67,74,86],[63,73,82],[59,72,78],
                                [76,91,95],[80,92,99],[84,93,103],[88,94,107],[97,110,115],[101,111,119],[105,112,123],
                                [109,113,127]])

    need_bridges=false

    generators = Vector{PauliBoolVec{N}}([])
    parameters = Vector{Float64}([])

    # Loop over trotter steps
    for ki in 1:k        
        # e^{i π/2 P2} e^{i π P1 /2}|ψ>
        ## ZZ layer
        for qubit in sequences
            for i in qubit[1]:qubit[2]
                pi = PauliBoolVec(N, Z=[i, i + 1])
                push!(generators, pi)
                push!(parameters, π)
            end
        end
        #bridges
        if need_bridges
            for link in bridges
                pi = PauliBoolVec(N, Z=[link[1], link[2]])
                push!(generators, pi)
                push!(parameters, π )
                pi = PauliBoolVec(N, Z=[link[2], link[3]])
                push!(generators, pi)
                push!(parameters, π )
            end
        else #PBC
            pi = PauliBoolVec(N,Z=[1])
            push!(generators, pi)
            push!(parameters, π )
        end   
        ## X layer
        # e^{i αn Pn / 2}
        for i in 1:N
            pi = PauliBoolVec(N, X=[i])
            push!(generators, pi)
            push!(parameters, α)
        end
    end
    return(generators,parameters)
end


function run(;α=.01, k=10, nsamples=10000, seed=1)
    
    N = 127  # number of qubits

    generators,parameters = eagle_processor(N,k,α)
#    print(generators)
    # o = PauliBoolVec(N,X=[1],Y=[2],Z=[3])

    #Mz
    # o = PauliBoolVec(N, X=[13,29,31], Y=[9,30], Z=[8,12,17,28,32])
    # o = PauliBoolVec(N, Z=[1,2,3,4])
    # o = PauliBoolVec(N, Y=[1], Z=[2,3,4])
    # o = PauliBoolVec(N, Z=[1,30])
    o = PauliBoolVec(N, Z=[1])


    ket = zeros(Bool, N)
    Random.seed!(seed)
    e = zeros(ComplexF64, nsamples)
    # e2 = zeros(ComplexF64, nsamples)
    
    ei, _ = UnitaryPruning.stochastic_pauli_dynamics_run(generators, parameters, o, ket)
    e[1] = ei 
    
    for i in 2:nsamples
        ei, _ = UnitaryPruning.stochastic_pauli_dynamics_run(generators, parameters, o, ket)
        e[i] = (e[i-1]*(i-1) + ei)/i
    end
    return e
end


function run_even(;α=.01, k=10, sample_size=10, nsamples=10000, seed=1)
    
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

    # our operator
    o = PauliBoolVec(N, Z=[1])


    ket = zeros(Bool, N)
    Random.seed!(seed)
 
    sample_set = zeros(nsamples)
    for i in 1:nsamples
        e = 0.0
        for _ in 1:sample_size
            e += UnitaryPruning.get_random_leaf(generators, parameters, o, ket)[1]
        end
        sample_set[i] = e
    end
    return mean(sample_set), std(sample_set) 
end


nruns = 100
nsamples = 10000

final_vals_stoc = []
final_vals_errs = []

for i in 1:16
    trajectories = []
    avg_traj = zeros(nsamples) 
    for runi in 1:nruns
        ei = run(α=i*π/32, k=5, nsamples=nsamples, seed=runi)
        ei = real(ei)
        push!(trajectories, ei)

        avg_traj .+= ei
    end
    avg_traj ./= nruns
    var_traj = sum([(i .- avg_traj).*(i .- avg_traj) for i in trajectories]) ./nruns
    
    std_traj = sqrt.(var_traj)
    z = 1.96 # 95% confidence
    std_traj = z .* std_traj ./ sqrt(nruns)

    @printf(" α: %4i avg: %12.8f ± %-12.6f var: %12.8f\n", i, avg_traj[end], std_traj[end], var_traj[end])

    plot(trajectories, ylim=[-1,1], legend=false, color=:grey, alpha=.1)
    plot!(avg_traj, color=:black, ribbon=std_traj)
    # plot!(avg_traj .+ std_traj, color=:red)
    # plot!(avg_traj .- std_traj, color=:red)
    savefig(@sprintf "%05i.pdf" i)
    push!(final_vals_stoc, avg_traj[end])
    push!(final_vals_errs, std_traj[end])
end

plot(final_vals_stoc, ribbon=final_vals_errs, marker=:circle); 
#plot!(final_vals); 
savefig("plot.pdf")
