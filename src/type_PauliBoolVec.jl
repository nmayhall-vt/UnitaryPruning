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
    θ += 2*count(p1.x .& p2.z)
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
    Base.display(p::PauliBoolVec{N}) where {N}

TBW
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
                str = str * "iY"
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