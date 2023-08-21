"""
In this representation, the Pauli string operator is represented as two binary strings, one for x and one for z.

The format is as follows: Z^z1 X^x1 ⊗ Z^z2 X^x2 ⊗ ⋯ ⊗ Z^zN X^xN  
    
Products of operators simply concatonate the left and right strings separately. For example, 

    XYZIY = 11001|01101
"""
struct PauliBitString{N} <: Pauli 
    θ::UInt8
    z::Int128
    x::Int128
end


"""
    PauliBitString(z::I, x::I) where I<:Integer

TBW
"""
function PauliBitString(z::I, x::I) where I<:Integer
    N = maximum(map(i -> ndigits(i, base=2), [x, z]))
    θ = count_ones(z & x)*3 % 4
    return PauliBitString{N}(θ, z, x)
end

function PauliBitString(z::I, x::I, N) where I<:Integer
    θ = count_ones(z & x)*3 % 4
    return PauliBitString{N}(θ, z, x)
end

"""
    PauliBitString(str::String)

TBW
"""
function PauliBitString(str::String)
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
    θ = 3*ny%4
    return PauliBitString{N}(θ, z,x) 
end


"""
    PauliBoolVec(N, X::Vector, Y::Vector, Z::Vector)

constructor for creating PauliBoolVec by specifying the qubits where each X, Y, and Z gates exist 
"""
function PauliBitString(N::Integer; X=[], Y=[], Z=[])
    for i in X
        i ∉ Y || throw(DimensionMismatch)
        i ∉ Z || throw(DimensionMismatch)
    end
    for i in Y
        i ∉ Z || throw(DimensionMismatch)
    end
    
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
    return PauliBitString(join(str))
    
end



"""
    multiply(p1::PauliBitString{N}, θ1,  p2::PauliBitString{N}, θ2) where N

TBW
"""
function multiply(p1::PauliBitString{N},  p2::PauliBitString{N}) where N
#    print("p1 ", p1, " p2 ", p2)
#    print("p1.x ", p1.x, "\n", "p2.x ", p2.x, "\n")
#    print("p1.z ", p1.z, "\n", "p2.z ", p2.z, "\n")
    x = p1.x ⊻ p2.x
    z = p1.z ⊻ p2.z
#    print("p1.θ ", p1.θ, "\n", "p2.θ ", p2.θ, "\n")
#    exit()
    θ = (p1.θ + p2.θ ) % 4
    θ = (θ + 2*count_ones(p1.x & p2.z)) % 4
    return PauliBitString{N}(θ,z,x)
end

"""
    Base.:*(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}

TBW
"""
function Base.:*(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}
    return multiply(p1,p2)
end


"""
    is_diagonal(p::PauliBitString)

Does this operator consist of only I and/or Z?
"""
function is_diagonal(p::PauliBitString)
    return count_ones(p.x) == 0
end

"""
    Base.show(io::IO, P::PauliMask)

TBW
"""
function Base.show(io::IO, p::PauliBitString{N}) where N
    # print(io, @sprintf "Pstring(P))
    # println(io, 1im^p.θ,"|", string(p)) 
    println(@sprintf "%2i %2iim | %s" real(1im^p.θ) imag(1im^p.θ) string(p)) 
end


function Base.:+(p::PauliBitString{N}, θ::Integer) where N
    return PauliBitString{N}((p.θ + θ)%4, p.z, p.x)
end

"""
    Base.display(p::PauliBitString)

Display, y = iY
"""
function Base.string(p::PauliBitString{N}) where N
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
    random_PauliBitString(N)

TBW
"""
function random_PauliBitString(N)
    return PauliBitString{N}(rand(0:3), rand(Int128),rand(Int128))
end

"""
    Base.:(==)(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}

Check if they are equal, return true or false
"""
function Base.:(==)(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}
    return p1.x == p2.x && p1.z == p2.z && p1.θ == p2.θ
end


"""
    commute(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}

Check if they commute, return true or false
"""
function commute(p1::PauliBitString{N}, p2::PauliBitString{N}) where {N}
    return iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 
end


"""
    expectation_value_sign(p::PauliBitString{N}, ket::Vector{Bool}) where N

compute expectation value of PauliBitString `o` for a product state `ket`
"""
function expectation_value_sign(p::PauliBitString{N}, ket::BasisState{N}) where N
    is_diagonal(p) || return 0.0

    count_ones(p.z & ket.v) % 2 == 0 || return -(1im)^p.θ
    return (1im)^p.θ 
end
