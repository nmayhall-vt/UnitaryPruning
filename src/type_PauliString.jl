using StaticArrays

struct PauliString{N}
    string::MVector{N,Char}
end

function PauliString(str::String) 
    return PauliString(MVector{length(str)}([i for i in str]))
end

function Base.print(ps::PauliString{N}) where N 
    for i in ps.string
        print(i)
    end
    println()
end
