import Base
#using LinearAlgebra
import IterTools: enumerate

########################
# Samantha Barren '22
# modified
########################

struct PauliMask{T<:Unsigned, N} <: Pauli
    #=
    Pauli mask representation of a Pauli string.
    This type is parametric on `T` so that smaller
    sized integers can be used to represent the Pauli
    string, if possible. This also performs a check that
    the given integers constitute a valid Pauli string.

    Example:
    Pauli string: "IIXYZYZIYIZ"
              id:  11000001010 -> 1546
               x:  00100000000 -> 256
               y:  00010100100 -> 164
               z:  00001010001 -> 81

    =#
    id::T
    x::T
    y::T
    z::T
end

function PauliMask(str::String, T=UInt64)
    return pauli_string_to_pauli(str, T)
end
function PauliMask(x::T, y::T, z::T, N)  where T
    return PauliMask{T,N}(~(x|y|z), x, y, z)
end


Base.:(==)(lhs ::PauliMask, rhs ::PauliMask) = (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z) && (lhs.id == rhs.id)


function num_qubits(P::PauliMask)
    return maximum(map(k -> ndigits(k, base=2), [P.x, P.y, P.z]))
end


function pauli_to_axes(P::PauliMask{T,N}) where {T,N}
    if num_qubits(P) > N
        error("Needs more qubits")
    end

    l = zeros(Int64, N)

    for i=1:N
        for (ax, p_comp) in zip([1, 2, 3], [P.x, P.y, P.z])
            res = ((p_comp >> (i-1)) & 1)
            if res == 1
                l[i] = ax
            end
        end
    end
    return reverse(l)
end


function pauli_string_to_pauli(ps::String; T = UInt64)
    l = zeros(Int64, length(ps))
    for (i, c)=enumerate(ps)
        if c == 'I'
            l[i] = 0
        elseif c == 'X'
            l[i] = 1
        elseif c == 'Y'
            l[i] = 2
        elseif c == 'Z'
            l[i] = 3
        else
            error("Invalid character: $c")
        end
    end
    return pauli_string_to_pauli(l, T=T)
end


function pauli_string_to_pauli(ps::Vector{<:Integer}; T = UInt64)
    pm = [0, 0, 0, 0]
    _pauli_masks(pm, ps)
    _type_out = unsigned(T)
    return PauliMask(
        _type_out(pm[2]),
        _type_out(pm[3]),
        _type_out(pm[4]),
        length(ps)
    )
end

function pauli_to_pauli_string(P::PauliMask{T,N}) where {T,N}
    plist = ["I","X","Y","Z"]
    pax = pauli_to_axes(P) .+ 1
    pstr = []
    for el in reverse(pax)
        push!(pstr,plist[el])
    end
    return join(pstr)
end

function _pauli_masks(res::Array{Int64,1}, pauli_str::Array{Int64,1})
    for (i,ax)=enumerate(reverse(pauli_str))
        res[ax+1] += 2^(i-1)
    end
end


function Base.show(io::IO, P::PauliMask)
    num_qubits = maximum(map(i -> ndigits(i, base=2), [P.x, P.y, P.z]))
    xs = bitstring(P.x)[end-num_qubits+1:end]
    ys = bitstring(P.y)[end-num_qubits+1:end]
    zs = bitstring(P.z)[end-num_qubits+1:end]
    print(io, "Pauli(x=$xs, y=$ys, z=$zs)")
end


#function phase_shift(alpha::ComplexF64, i::Integer)
#    if i == 0
#        return     alpha.re + im*alpha.im
#    elseif i == 1
#        return  im*alpha.re - alpha.im
#    elseif i == 2
#        return    -alpha.re - im*alpha.im
#    else
#        return -im*alpha.re + alpha.im
#    end
#end


#function pauli_phase(pm::PauliMask, a::Int64)
#    # Compute the phase gamma where P|a> = gamma |b> for
#    # a pauli string P and basis state |a>.
#    # Convention:
#    #   0 -> 1
#    #   1 -> +i
#    #   2 -> -1
#    #   3 -> -i
#    # pm = pauli_mask input
#    x = count_ones((pm.y | pm.z) & a) % 2
#    y = count_ones(pm.y) % 4
#
#    alpha = y
#    beta = 2*x
#
#    return UInt8((alpha+beta) % 4)
#end


#function pauli_apply(pm::PauliMask, a::Int64)
#    # pm = pauli_mask input
#    return xor((pm.x|pm.y),a)
#end


#function pauli_mult!(pm::PauliMask, state::Array{ComplexF64,1}, result::Array{ComplexF64,1})
#    N = length(state)
#    for i=0:N-1
#        j = pauli_apply(pm, i)
#        phase = pauli_phase(pm, i) % 4
#        r = state[i+1]
#        result[j+1] = phase_shift(r, phase)
#    end
#end


function commute(P::PauliMask, Q::PauliMask)
    id = (P.id | Q.id)
    x = (P.x & Q.x)|id
    y = (P.y & Q.y)|id
    z = (P.z & Q.z)|id

    return iseven(count_ones(x | y | z))

    #res = 0
    #res += count_ones(x)
    #res += count_ones(y)
    #res += count_ones(z)
    #return Bool((res+1)%2)
end


function commutator(P::PauliMask{T,N}, Q::PauliMask{T,N}) where {T,N}
    phase = T(0)

    x = (P.z & Q.y) | (Q.z & P.y) | (P.id & Q.x) | (Q.id & P.x)
    phase += 3*(count_ones(P.z & Q.y) - count_ones(Q.z & P.y))

    y = (P.x & Q.z) | (Q.x & P.z) | (P.id & Q.y) | (Q.id & P.y)
    phase -= 1*(count_ones(P.x & Q.z) - count_ones(Q.x & P.z))

    z = (P.x & Q.y) | (Q.x & P.y) | (P.id & Q.z) | (Q.id & P.z)
    phase += 1*(count_ones(P.x & Q.y) - count_ones(Q.x & P.y))

    phase = phase%4
    if phase == 0
        return 1,  PauliMask{T,N}(~(x|y|z), x, y, z)
    elseif phase == 1
        return 1im,  PauliMask{T,N}(~(x|y|z), x, y, z)
    elseif phase == 2
        return -1,  PauliMask{T,N}(~(x|y|z), x, y, z)
    elseif phase == 3
        return -1im,  PauliMask{T,N}(~(x|y|z), x, y, z)
    else
        throw(DomainError)
    end
    #return 1im^phase, PauliMask{T,N}(~(x|y|z), x, y, z)
end


"""
    expectation_value_sign(o::PauliMask{T,N}, ket::Vector{Integer})

compute expectation value of PauliString `o` for a product state `ket`
"""
function expectation_value_sign(pm::PauliMask{T,N}, ket::T) where {T,N}
    #is_diagonal(pm) || return 0.0
    #return (-1)^count_ones(pm.z & ket)%2
    return iseven(count_ones(pm.z & ket))
end


"""
    is_diagonal(pm1::PauliMask{T,N})

Is `pm1` diagonal?
"""
function is_diagonal(pm1::PauliMask{T,N}) where {T,N}
    return count_ones(pm1.x | pm1.y) == 0 
end

