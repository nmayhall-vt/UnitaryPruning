using StaticArrays
using Parameters

#
#   (i)^θ ⋅ Z^z1 X^x1 ⊗ Z^z2 X^x2 ⊗ ⋯ ⊗ Z^zN X^xN  
mutable struct PauliBoolVec{N} <: Pauli
    θ::UInt8
    x::MVector{N, Bool}
    z::MVector{N, Bool}
end

function PauliBoolVec(str::String)
    for i in str
        i in ['I', 'Z', 'X', 'Y'] || error("Bad string: ", str)
    end

    N = length(str)
    x = MVector{N, Bool}([i in ['X', 'Y'] for i in str])
    z = MVector{N, Bool}([i in ['Z', 'Y'] for i in str])
    θ = 3*count([i == 'Y' for i in str]) % 4
    return PauliBoolVec{N}(θ, x, z)
end


"""
    PauliBoolVec(N, X::Vector, Y::Vector, Z::Vector)

constructor for creating PauliBoolVec by specifying the qubits where each X, Y, and Z gates exist 
"""
function PauliBoolVec(N::Integer; X=[], Y=[], Z=[])
    for i in X
        i ∉ Y || throw(DimensionMismatch)
        i ∉ Z || throw(DimensionMismatch)
    end
    for i in Y
        i ∉ Z || throw(DimensionMismatch)
    end
    
    θ = 3*length(Y) % 4
    p = PauliBoolVec{N}(θ, zeros(N), zeros(N))
   
    for i in X
        p.x[i] = true
    end
    for i in Y
        p.x[i] = true
        p.z[i] = true
    end
    for i in Z
        p.z[i] = true
    end
    return p
end

function nick(;a=[])
    display(a)
end

"""
    multiply(p1::PauliBoolVec{N}, p2::PauliBoolVec{N}) where {N}

TBW
"""
function multiply(p1::PauliBoolVec{N}, p2::PauliBoolVec{N}) where {N}

    # x = x1 + x2 % 2
    # z = z1 + z2 % 2
    # sign = x2 + z2

    x = p1.x .⊻ p2.x
    z = p1.z .⊻ p2.z
    θ = (p1.θ + p2.θ ) % 4
    # println(p1.θ, " ", p2.θ)
    θ = (θ + 2*count(p1.x .& p2.z)) % 4
    # println(θ)
    return PauliBoolVec{N}(θ,x,z)
end

"""
    Base.:*(p1::PauliBoolVec{N}, p2::PauliBoolVec{N}) where {N}

TBW
"""
function Base.:*(p1::PauliBoolVec{N}, p2::PauliBoolVec{N}) where {N}
    return multiply(p1,p2)
end


"""
    Base.:(==)(p1::PauliBoolVec{N}, p2::PauliBoolVec{N}) where {N}

TBW
"""
function Base.:(==)(p1::PauliBoolVec{N}, p2::PauliBoolVec{N}) where {N}
    return (p1.θ == p2.θ) & all(p1.x .== p2.x) & all(p1.z .== p2.z) 
end


"""
    commute_check(p1::PauliBoolVec{N}, p2::PauliBoolVec{N}) where {N}

Do Pauli strings `p1` and `p2` commute?
"""
function commute(p1::PauliBoolVec{N}, p2::PauliBoolVec{N}) where {N}
    return iseven(count(p1.x .& p2.z) - count(p1.z .& p2.x)) 
end


function to_matrix(p::PauliBoolVec{N}, T=ComplexF64) where N
    x = [0 1; 1 0]
    y = [0 1; -1 0]
    z = [1 0; 0 -1]
    I = [1 0; 0 1]

    out = zeros(T, 2, 2)
    if p.x[1] == 1 && p.z[1] == 1
        out .= y
    elseif p.x[1] == 1 && p.z[1] == 0
        out .= x
    elseif p.x[1] == 0 && p.z[1] == 1
        out .= z
    elseif p.x[1] == 0 && p.z[1] == 0
        out .= I
    end
    
    for i in 2:N
        if p.x[i] == 1 && p.z[i] == 1 
            out = kron(out, y)
        elseif p.x[i] == 1 && p.z[i] == 0 
            out = kron(out, x)
        elseif p.x[i] == 0 && p.z[i] == 1
            out = kron(out, z)
        elseif p.x[i] == 0 && p.z[i] == 0
            out = kron(out, I)
        end
    end
    return (1im)^p.θ .* out
end


"""
    is_diagonal(p::PauliBoolVec{N}) where N

is `p` diagonal? (i.e., comprised of only I and Z operators)
"""
function is_diagonal(p::PauliBoolVec{N}) where N
    return count(p.x) == 0
end


"""
    Base.display(p::PauliBoolVec{N}) where {N}

Display, y = iY
"""
function Base.display(p::PauliBoolVec{N}) where {N}
    str = ""
    for i in 1:N
        if p.x[i] == 0 
            if p.z[i] == 0
                str = str * "I"
            elseif p.z[i] == 1
                str = str * "Z"
            end
        end 
        if p.x[i] == 1 
            if p.z[i] == 0
                str = str * "X"
            elseif p.z[i] == 1
                str = str * "y"
            end
        end
        if i < N
            str = str * "⊗"
        end 
    end
    println(1im^p.θ,"|", str) 
end

"""
    countY(p::PauliBoolVec{N}) where {N}

Return the number of Y pauli's in `p`
"""
function countY(p::PauliBoolVec{N}) where {N}
    return count(p.x .& p.z)
end

"""
    countI(p::PauliBoolVec{N}) where {N}

Return the number of I pauli's in `p`
"""
function countI(p::PauliBoolVec{N}) where {N}
    return count(p.x .⊽ p.z)
end

"""
    countX(p::PauliBoolVec{N}) where {N}

Return the number of X pauli's in `p`
"""
function countX(p::PauliBoolVec{N}) where {N}
    return count(p.x .& (p.x .⊻ p.z))
end

"""
    countZ(p::PauliBoolVec{N}) where {N}

Return the number of Z pauli's in `p`
"""
function countZ(p::PauliBoolVec{N}) where {N}
    return count(p.z .& (p.x .⊻ p.z))
end


"""
    expectation_value_sign(o::PauliString{N}, ket::Vector) where N

compute expectation value of PauliString `o` for a product state `ket`
"""
function expectation_value_sign(p::PauliBoolVec{N}, ket::Vector{Bool}) where N
    length(ket) == N || error(" ket and paulistring don't match") 
   
    is_diagonal(p) || return 0.0

    count(p.z .& ket) % 2 == 0 || return -(1im)^p.θ
    return (1im)^p.θ 
    
    # sign = 1
    # for i in 1:N
    #     if o.z[i] == 1
    #         if ket[i] == 1
    #             sign = -sign
    #         end
    #     end
    # end
    # return sign
end

"""
    build_time_evolution_matrix(gs::Vector{PauliBoolVec{N}}, angles::Vector)


Build the Unitary matrix: exp(i αn Pn / 2) ⋯ exp(i α2 P2 / 2) exp(i α1 P1 / 2)
"""
function build_time_evolution_matrix(generators::Vector{PauliBoolVec{N}}, angles::Vector) where {N}
    U = to_matrix(PauliBoolVec(N))
    Nt = length(generators)
    length(angles) == Nt || throw(DimensionMismatch)
    for t in 1:Nt
        α = angles[t]
        # Ut = e(i α Pn) = cos(α) I + i sin(α) Pn
        U = cos(α/2) .* U   .+   1im * sin(α/2) .* to_matrix(generators[t]) * U 
    end
    return U
end

