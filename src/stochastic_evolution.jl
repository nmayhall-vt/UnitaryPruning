using Random

"""
    stochastic_pauli_dynamics_run(generators::Vector{PauliBoolVec{N}}, angles, o::PauliBoolVec{N}) where N

TBW
"""
function stochastic_pauli_dynamics_run(generators::Vector{PauliBoolVec{N}}, angles, o_in::PauliBoolVec{N}, ket) where N

    o = deepcopy(o_in)
    #
    # for a single pauli Unitary, Un = exp(-i θn Pn/2)
    # U' O U = cos(θ/2) O - i sin(θ/2) PO
    nt = length(generators)
    length(angles) == nt || throw(DimensionMismatch)

    scale = 1.0
   
    path = zeros(Int,nt)
    # sin_angles = sin.(angles)
    # cos_angles = cos.(angles)
    scales = sin.(angles) .+ cos.(angles) 
    bias = sin.(angles) ./ scales
    # bias = sin.(angles).^2 
    # bias = tan.(angles) 
  
    for t in reverse(1:nt)
        g = generators[t]
        commute(o,g) == false || continue

        path[t] = -1
        
        if rand() < bias[t]
            # sin branch
           
            # println("sin: ", bias[t])
            # error("here")
            o = multiply(g,o)
            o.θ = (o.θ + 3) % 4

            path[t] = 1
        end
        scale *= scales[t]
        # scale *= cos_angles[t] + sin_angles[t]
        # display(o)
        # @printf("%4i ", t)
        # println(o) 
    end

    return scale * expectation_value_sign(o,ket), path

end