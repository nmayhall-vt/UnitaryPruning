using Random

"""
    stochastic_pauli_dynamics_run(generators::Vector{PauliBoolVec{N}}, angles, o::PauliBoolVec{N}) where N

TBW
"""
function stochastic_pauli_dynamics_run(generators::Vector{PauliBoolVec{N}}, angles, o_in::PauliBoolVec{N}, ket) where N

    o = deepcopy(o_in)
    #
    # for a single pauli Unitary, Un = exp(-i θn Pn/2)
    # U' O U = cos(θ/2) O + i sin(θ/2) PO
    nt = length(generators)
    length(angles) == nt || throw(DimensionMismatch)
    
    sin_params = sin.(angles./2)
    cos_params = cos.(angles./2)
    tan_params = tan.(angles./2) 
    # tan_params = tan.(angles.*2) .* tan.(angles.*2)
    # bias = tan.(parameters .+ 2\pi)
  
    for t in reverse(1:nt)
        g = generators[t]
        commute(o,g) == false || continue

        branch = 0 #cosine branch
        if rand() < tan_params[t]
            # sin branch
            o = multiply(o,g)
            o.θ = (o.θ + 1) % 2
        end
        # @printf("%4i ", t)
        # println(o) 
    end

    return expectation_value_sign(o,ket) * (1im)^o.θ

end