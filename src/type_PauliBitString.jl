#struct PauliBitString{N}
#    #x::MVector{N,Bool}
#    #z::MVector{N,Bool}
#    x::BitVector{N}
#    y::BitVector{N}
#end
#
#function PauliBitString(x::T, z::T, N::T) where {T<:Integer}
#    return PauliBitString{N}(digits(x, base=2, pad=N), digits(z, base=2, pad=N))
#end
#
#function PauliBitString(ps::PauliString{N}) where {N}
#    #x = [0 for i in 1:N]
#    #z = [0 for i in 1:N]
#    #println(join((ps[i] == 'X' || ps[i] == 'Y' ? true : false  for i in 1:N)))
#    #println(join((ps[i] == 'Z' || ps[i] == 'Y' ? true : false  for i in 1:N)))
#    #println(((ps[i] == 'I' || ps[i] == 'Z') ? 1 : for i in ps.string))
#    return PauliBitString(MVector{N,Bool}(ps[i] == 'X' || ps[i] == 'Y' ? true : false  for i in 1:N),
#                          MVector{N,Bool}(ps[i] == 'Z' || ps[i] == 'Y' ? true : false  for i in 1:N))
#end
#
##function test4(a, b)
##    s = length(a)
##    c = Vector{Bool}(undef, s)
##    @avx for i in 1:s
##        c[i] = a[i] > b[i]
##    end
##    c
##end



"""
    PauliBitString{N} <: Operator 

In this representation, the Pauli string operator is represented as a binary string.
Any possible Pauli string can be represented in this way, referred to as the symplectic representation
which finds use in quantum error correction. 

    
    I : 0|0 
    Z : 0|1 
    X : 1|0 
    Y : 1|1 

Products of operators simply concatonate the left and right strings separately. For example, 

    XYZIY = 11001|01101

"""
struct PauliBitString{T, N} 
    x::T  
    z::T
    #function PauliBitString(str::String) where {N}
    #    value = parse(U, string("0b", str))
    #    nbits = length(str)
    #    return new{nbits}(value)
    #end
end


#Base.print(pbs::PauliBitString) = print(pbs.x, pbs.z)
Base.println(pbs::PauliBitString) = println(string(pbs))
Base.string(pbs::PauliBitString) = join([bitstring(pbs.x), "|", bitstring(pbs.z)]) 

function Base.print(pbs::PauliBitString{T,N}) where {T,N}
    for i in N:1
        print(i)
    end
end


"""
    is_diagonal(pbs::PauliBitString)

Does this operator consist of only I and/or Z?
"""
function is_diagonal(pbs::PauliBitString)
    return pbs.x == 0
end

"""
    n_non_eye(pbs::PauliBitString{T,N}) where {T,N}

Return number of non-identity (i.e., X, Y, or Z) in string
"""
function n_non_eye(pbs::PauliBitString{T,N}) where {T,N}
    return count_ones(pbs.x | pbs.z)
end

"""
"""
#@inline function commute(pbs1::PauliBitString{T,N}, pbs2::PauliBitString{T,N}) where {T,N}
#function commute(pbs1::PauliBitString{T,N}, pbs2::PauliBitString{T,N}) where {T,N}
function commute(pbs1, pbs2)
    return iseven(count_ones(pbs1.x & pbs2.z) - count_ones(pbs1.z & pbs2.x)) 
    #return iseven(sum(pbs1.x & pbs2.z) - sum(pbs1.z & pbs2.x)) 
end


"""
    commutator(pbs1::PauliBitString, pbs2::PauliBitString)
"""
function commutator(pbs1::PauliBitString{T,N}, pbs2::PauliBitString{T,N}) where {T,N}
    #return 1im^(count_ones(pbs1.x & pbs2.z)),  PauliBitString{T,N}(pbs1.x ⊻ pbs2.x, pbs1.z ⊻ pbs2.z) 
    #return -1im*1im^(count_ones(pbs1.x & pbs2.z) - count_ones(pbs1.z & pbs2.x)) * 1im^(count_ones(pbs1.x & pbs2.z)),  PauliBitString{T,N}(pbs1.x ⊻ pbs2.x, pbs1.z ⊻ pbs2.z) 
    #return 1im^(-count_ones(pbs1.z & pbs2.x)),  PauliBitString{T,N}(pbs1.x ⊻ pbs2.x, pbs1.z ⊻ pbs2.z) 
end


function string_to_bits(str::String)
    x = zeros(UInt8, length(str))
    z = zeros(UInt8, length(str))
    for (i_idx,i) in enumerate(str)
        if i == "X"
            x[i_idx] = 1
        elseif i == "Y" 
            x[i_idx] = 1
            z[i_idx] = 1
        elseif i == "Z" 
            z[i_idx] = 1
        end
    end
    return x,z
end


