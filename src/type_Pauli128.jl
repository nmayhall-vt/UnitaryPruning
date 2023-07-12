"""
    Pauli128 <: Operator 

In this representation, the Pauli string operator is represented as two binary strings, one for x and one for z.
Here, the phase is not part of the type, that way we can use it as a Dict key. This means that some functions must behave differently.
E.g., multiply, must return not only a new `Pauli64`, but also the associated phase the results.

The format is as follows: Z^z1 X^x1 ⊗ Z^z2 X^x2 ⊗ ⋯ ⊗ Z^zN X^xN  
    
Products of operators simply concatonate the left and right strings separately. For example, 

    XYZIY = 11001|01101

"""
struct Pauli128{N} <: Pauli 
    z::Int128
    x::Int128
end


"""
    Pauli128(z::I, x::I) where I<:Integer

TBW
"""
function Pauli128(z::I, x::I) where I<:Integer
    N = maximum(map(i -> ndigits(i, base=2), [x, z]))
    return Pauli128{N}(z,x), 1
end

function Pauli128(str::String)
    for i in str
        i in ['I', 'Z', 'X', 'Y'] || error("Bad string: ", str)
    end

    x = Int128(0)
    z = Int128(0)
    ny = 0 
    N = length(str)
    idx = Int128(1)
    two = Int128(2)
    one = Int128(1)

    for i in str
        # println(i, " ", idx, typeof(idx))
        if i in ['X', 'Y']
            x |= two^(idx-one)
            if i == 'Y'
                ny += 1
            end
        end
        if i in ['Z', 'Y']
            z |= two^(idx-one)
        end
        idx += 1
    end
    return Pauli128{N}(z,x), 3*ny%4
end


"""
    PauliBoolVec(N, X::Vector, Y::Vector, Z::Vector)

constructor for creating PauliBoolVec by specifying the qubits where each X, Y, and Z gates exist 
"""
function Pauli128(N::Integer; X=[], Y=[], Z=[])
    for i in X
        i ∉ Y || throw(DimensionMismatch)
        i ∉ Z || throw(DimensionMismatch)
    end
    for i in Y
        i ∉ Z || throw(DimensionMismatch)
    end
    
    θ = 3*length(Y) % 4

    str = ["I" for i in 1:N]
    for i in X
        str[i] = "X"
    end
    for i in Y
        str[i] = "Y"
    end
    for i in Z
        str[i] = "Z"
    end
   
    # print(str[1:N])
    return Pauli128(join(str))
    
end


##############################################################################

struct Ket128
    v::Int128
end
           
function Ket128(vec::Vector)
    one = Int128(1)
    two = Int128(2)

    v = Int128(0)
    for i in 1:length(vec)
        if vec[i] == 1
            v |= two^(i-one)
        end
    end
    return Ket128(v)
end
##############################################################################

"""
    multiply(p1::Pauli128{N}, θ1,  p2::Pauli128{N}, θ2) where N

TBW
"""
function multiply(p1::Pauli128{N}, θ1,  p2::Pauli128{N}, θ2) where N
    x = p1.x ⊻ p2.x
    z = p1.z ⊻ p2.z
    θ = (θ1 + θ2 ) % 4
    θ = (θ + 2*count_ones(p1.x & p2.z)) % 4
    return Pauli128{N}(z,x), θ
end

"""
    is_diagonal(p::Pauli128)

Does this operator consist of only I and/or Z?
"""
function is_diagonal(p::Pauli128)
    return count_ones(p.x) == 0
end

"""
    Base.show(io::IO, P::PauliMask)

TBW
"""
function Base.show(io::IO, P::Pauli128{N}) where N
    # xs = bitstring(P.x)[end-N+1:end]
    # zs = bitstring(P.z)[end-N+1:end]
    # print(io, "Pauli(z=$zs, x=$xs)")
    print(io, string(P))
end


"""
    Base.display(p::Pauli128)

Display, y = iY
"""
function Base.string(p::Pauli128{N}) where N
    Iloc = get_on_bits(p.x ⊽ p.z)
    yloc = get_on_bits(p.x & p.z)
    Xloc = get_on_bits(p.x & ~p.z)
    Zloc = get_on_bits(p.z & ~p.x)
    out = ["I" for i in 1:128]

    for i in Xloc
        out[i] = "X"
    end
    for i in yloc
        out[i] = "y"
    end
    for i in Zloc
        out[i] = "Z"
    end
    return join(out[1:N])
end


"""
    random_Pauli128(N)

TBW
"""
function random_Pauli128(N)
    return Pauli128{N}(rand(Int128),rand(Int128))
end

"""
    Base.:(==)(p1::Pauli128{N}, p2::Pauli128{N}) where {N}

Check if they are equal, return true or false
"""
function Base.:(==)(p1::Pauli128{N}, p2::Pauli128{N}) where {N}
    return p1.x == p2.x && p1.z == p2.z
end


"""
    commute(p1::Pauli128{N}, p2::Pauli128{N}) where {N}

Check if they commute, return true or false
"""
function commute(p1::Pauli128{N}, p2::Pauli128{N}) where {N}
    return iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 
end


"""
    expectation_value_sign(p::Pauli128{N}, ket::Vector{Bool}) where N

compute expectation value of Pauli128 `o` for a product state `ket`
"""
function expectation_value_sign(p::Pauli128{N}, ket::Int128) where N
    is_diagonal(p) || return 0.0

    count_ones(p.z & ket) % 2 == 0 || return -1
    return 1 
end

"""
    expectation_value_sign(p::Tuple{Pauli128{N}, Int}, ket::Int128) where N

compute expectation value of Pauli128 `o` for a product state `ket`
    `p[1]` is the `Pauli128` and `p[2]` is the associated phase
"""
function expectation_value_sign(p::Tuple{Pauli128{N}, Int}, ket::Int128) where N
    is_diagonal(p[1]) || return 0.0

    count_ones(p[1].z & ket) % 2 == 0 || return -(1im)^p[2]
    return (1im)^p[2] 
end