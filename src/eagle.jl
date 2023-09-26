

function eagle_processor(o::Pauli{N}; k=10, α=0.01, sequences=[[1, 13], [19, 32], [38, 51], [57, 70], [76, 89], [95, 108], [114, 126]],
    bridges=[[1, 15, 19], [5, 16, 23], [9, 17, 27], [13, 18, 31], [33, 37, 52], [29, 36, 48], [25, 35, 44], [21, 34, 40],
        [38, 53, 57], [42, 54, 61], [46, 55, 65], [50, 56, 69], [71, 75, 90], [67, 74, 86], [63, 73, 82], [59, 72, 78],
        [76, 91, 95], [80, 92, 99], [84, 93, 103], [88, 94, 107], [97, 110, 115], [101, 111, 119], [105, 112, 123],
        [109, 113, 127]]) where {N}

    N == 127 || throw(DimensionMismatch)
    need_bridges = true

    generators = Vector{Pauli{N}}([])
    parameters = Vector{Float64}([])

    # Loop over trotter steps
    for ki in 1:k
        # e^{i π/2 P2} e^{i π P1 /2}|ψ>
        ## ZZ layer
        for qubit in sequences
            for i in qubit[1]:qubit[2]
                pi = Pauli(N, Z=[i, i + 1])
                push!(generators, pi)
                push!(parameters, π/2)
            end
        end
        #bridges
        if need_bridges
            for link in bridges
                pi = Pauli(N, Z=[link[1], link[2]])
                push!(generators, pi)
                push!(parameters, π/2)

                pi = Pauli(N, Z=[link[2], link[3]])
                push!(generators, pi)
                push!(parameters, π/2)
            end
        else #PBC
            pi = Pauli(N, Z=[1])
            push!(generators, pi)
            push!(parameters, π/2)
        end
        ## X layer
        # e^{i αn Pn / 2}
        for i in 1:N
            pi = Pauli(N, X=[i])
            pi = Pauli{N}((pi.θ + 2)%4, pi.pauli) # this accounts for the fact that the papers have -X and positive ZZ
            push!(generators, pi)
            push!(parameters, α)
        end
    end
    return (generators, parameters)
end
