"""
An occupation number vector, up to 128 qubits
"""
struct BasisState{N} 
    v::Int128
end
           
"""
    BasisState(vec::Vector{T}) where T<:Union{Bool, Integer}

TBW
"""
function BasisState(vec::Vector{T}) where T<:Union{Bool, Integer}
    two = Int128(2)
    v = Int128(0)

    for i in 1:length(vec)
        if vec[i] == 1
            v |= two^(i-1)
        end
    end
    return BasisState{N}(v)
end

"""
    BasisState(N::Integer, v::Integer)

TBW
"""
function BasisState(N::Integer, v::Integer)
    for i in N:128
        v &= ~(Int128(2)^(i-1))    
    end
    return BasisState{N}(Int128(v))
end


"""
    Base.show(io::IO, P::PauliBitString{N}) where N

TBW
"""
function Base.show(io::IO, P::BasisState{N}) where N
    print(io, string(P))
end

"""
    Base.display(p::Pauli128)

Display, y = iY
"""
function Base.string(p::BasisState{N}) where N
    out = [0 for i in 1:128]
    for i in get_on_bits(p.v)
        out[i] = 1
    end
    return join(out[1:N])
end
