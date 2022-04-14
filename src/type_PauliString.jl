using StaticArrays

struct PauliString{N} <: Pauli
    string::MVector{N,Char}
end

function PauliString(str::String)
    for i in str
        i in ['I', 'Z', 'X', 'Y'] || error("Bad string: ", str)
    end

    return PauliString(MVector{length(str)}([i for i in str]))
end

function PauliString(n::Int)
    return PauliString(MVector{n}(('I' for i in 1:n)))
end


function Base.println(ps::PauliString{N}, h::Real) where N 
    print(ps)
    @printf(" %14.10f\n", h)
end

Base.println(ps::PauliString{N}) where N = println(join(ps.string)) 
Base.print(ps::PauliString{N}) where N = print(join(ps.string)) 
Base.string(ps::PauliString) = join(ps.string) 
Base.length(ps::PauliString) = length(ps.string)
Base.getindex(ps::PauliString, i) = ps.string[i] 
Base.setindex!(ps::PauliString,val,key) = ps.string[key] = val
Base.:(==)(x::PauliString{N}, y::PauliString{N}) where N = all(x.string .== y.string);

function commute(ps1::PauliString{N}, ps2::PauliString{N}) where N
    #iseven(sum(ps1.string .== ps2.string))
    n_noncomm::Int = 0
    for i in 1:N
        if ps1[i] == 'I' || ps2[i] == 'I' 
            continue 
        end
        if ps1[i] == ps2[i]
            continue
        end
        n_noncomm += 1
    end
    return iseven(n_noncomm)
end

"""
    commutator(ps1::PauliString{N}, ps2::PauliString{N})
"""
function commutator(ps1::PauliString{N}, ps2::PauliString{N}) where N
#={{{=#
    # the number of non commuting terms gives the number of imaginary numbers
    # the product of the orderings gives the sign, e.g. [X..., Z...] adds factor of -1
    #commute(ps1, ps2) == false || error(" don't compute commutator for things that commute") 
    #commute == false || return 0, PauliString(0)
    ps3 = PauliString(N)
    phase::Complex{Int64} = 1
    #phase = 1
    for i::Int in 1:N
        if ps1[i] == 'I'
            ps3[i] = ps2[i]
        elseif ps1[i] == 'X'
            if ps2[i] == 'I'
                ps3[i] = 'X'
            elseif ps2[i] == 'X'
                ps3[i] = 'I'
            elseif ps2[i] == 'Z'
                phase *= -1im
                ps3[i] = 'Y'
            elseif ps2[i] == 'Y'
                phase *= 1im
                ps3[i] = 'Z'
            end
        elseif ps1[i] == 'Y'
            if ps2[i] == 'I'
                ps3[i] = 'Y'
            elseif ps2[i] == 'X'
                ps3[i] = 'Z'
                phase *= -1im
            elseif ps2[i] == 'Y'
                ps3[i] = 'I'
            elseif ps2[i] == 'Z'
                ps3[i] = 'X'
                phase *= 1im
            end
        elseif ps1[i] == 'Z'
            if ps2[i] == 'I'
                ps3[i] = 'Z'
            elseif ps2[i] == 'X'
                ps3[i] = 'Y'
                phase *= 1im
            elseif ps2[i] == 'Y'
                ps3[i] = 'X'
                phase *= -1im
            elseif ps2[i] == 'Z'
                ps3[i] = 'I'
            end
        end
    end
    return phase, ps3
    #return phase*2, ps3
end#=}}}=#



"""
    is_diagonal(ps1::PauliString{N})

Is `ps1` diagonal?
"""
function is_diagonal(ps1::PauliString{N}) where N
    for i in 1:N
        if ps1[i] == 'X' || ps1[i] == 'Y'
            return false
        end
    end
    return true
end


"""
    expectation_value_sign(o::PauliString{N}, ket::Vector{Integer})

compute expectation value of PauliString `o` for a product state `ket`
"""
function expectation_value_sign(o::PauliString{N}, ket::Vector) where N
    length(ket) == N || error(" ket and paulistring don't match") 
   
    is_diagonal(o) || return 0.0

    sign = 1
    for i in 1:N
        if o[i] == 'Z' && ket[i] == 1
            sign = -sign
        end
    end

    return sign
end


