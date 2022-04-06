using StaticArrays

struct PauliString{N}
    string::MVector{N,Char}
end

function PauliString(str::String)
    for i in str
        i in ['I', 'Z', 'X', 'Y'] || error("Bad string: ", str)
    end

    return PauliString(MVector{length(str)}([i for i in str]))
end

function Base.print(ps::PauliString{N}) where N 
    for i in ps.string
        print(i)
    end
    println()
end

Base.getindex(ps::PauliString, i) = ps.string[i] 
Base.setindex!(ps::PauliString,val,key) = ps.string[key] = val

function commute(ps1::PauliString{N}, ps2::PauliString{N}) where N
    #iseven(sum(ps1.string .== ps2.string))
    n_noncomm = 0
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
