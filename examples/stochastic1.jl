using UnitaryPruning
using Plots

function run(;α=.1, k=2)
    N = 6
    generators = Vector{PauliBoolVec{N}}([])
    parameters = Vector{Float64}([])

    # Loop over trotter steps
    for ki in 1:k
        # e^{i π/2 P2} e^{i π/2 P1}|ψ>
        ## ZZ layer
        for i in 1:N-1
            pi = PauliBoolVec(N, Z=[i, i + 1])
            push!(generators, pi)
            push!(parameters, π / 2)
        end

        ## X layer
        # e^{i αn Pn}
        for i in 1:N
            pi = PauliBoolVec(N, X=[i])
            push!(generators, pi)
            push!(parameters, α)
        end
    end

    # display.(generators)
    # display.(parameters)

    # o = PauliBoolVec(N,X=[1],Y=[2],Z=[3])

    #Mz
    o = to_matrix(PauliBoolVec(N, Z=[1]))
    for i in 2:N
        o .+= to_matrix(PauliBoolVec(N, Z=[i]))
    end

    tan_params = tan.(parameters)
    sin_params = sin.(parameters)
    cos_params = cos.(parameters)

    U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)

    m = diag(U*o*U')
    # println(" expectation values:")
    # display(m)
    return real(m[1]) 
end

vals = [];
for i in 1:20
    push!(vals, run(α=0.1, k=i))
end

plot(vals, marker = :circle)