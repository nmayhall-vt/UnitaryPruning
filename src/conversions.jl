"""
    PauliBitString(ps::PauliString{N}; T=UInt) where N
"""
function PauliBitString(ps::PauliString{N}; T=UInt) where N
    x = join(collect(ps[i] == 'X' || ps[i] == 'Y' ? 1 : 0 for i in reverse(1:N)))
    z = join(collect(ps[i] == 'Z' || ps[i] == 'Y' ? 1 : 0 for i in reverse(1:N)))
    return PauliBitString{T,N}(parse(T, x; base=2), parse(T, z; base=2))
end



"""
    PauliString(pbs::PauliBitString{T,N}) where {T,N}
"""
function PauliString(pbs::PauliBitString{T,N}) where {T,N}
    str = []
    x = bitstring(pbs.x)
    z = bitstring(pbs.z)
    for i in 0:N-1
        if x[end-i] == '0' && z[end-i] == '0'
            push!(str, 'I')
        elseif x[end-i] == '0' && z[end-i] == '1'
            push!(str, 'Z')
        elseif x[end-i] == '1' && z[end-i] == '0'
            push!(str, 'X')
        elseif x[end-i] == '1' && z[end-i] == '1'
            push!(str, 'Y')
        else
            error(" Something wrong here")
        end
    end
    return PauliString(MVector{N}(str))
end
