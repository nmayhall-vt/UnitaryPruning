using UnitaryPruning, PyCall

np = pyimport("numpy")

ham_ops = Vector{PauliString{N}}()
ham_par = Vector{Float64}()
ansatz_ops = Vector{PauliString{N}}()
ansatz_par = Vector{Float64}()

Na_ops = Vector{PauliString{N}}()
Na_par = Vector{Float64}()
Nb_ops = Vector{PauliString{N}}()
Nb_par = Vector{Float64}()

for i in np.load(dir*"ansatz_par.npy")
    # Ansatz parametersshould be real 
    # openfermion tends to add phase to coefficients
    push!(ansatz_par, real(-1im * i))
end
for i in np.load(dir*"ham_par.npy")
    push!(ham_par, i)
end
for i in np.load(dir*"Na_par.npy")
    push!(Na_par, i)
end
for i in np.load(dir*"Nb_par.npy")
    push!(Nb_par, i)
end


do_mask = true 

if do_mask
    ansatz_ops = Vector{PauliMask{T,N}}()
    ham_ops    = Vector{PauliMask{T,N}}()
    Na_ops     = Vector{PauliMask{T,N}}()
    Nb_ops     = Vector{PauliMask{T,N}}()
    ref_state = parse(T, join(ref_state); base=2)
    for i in np.load("ham_ops.npy")
        push!(ham_ops, PauliMask(PauliString(i)))
    end
    for i in np.load(dir*"ansatz_ops.npy")
        push!(ansatz_ops, PauliMask(PauliString(i)))
    end
    for i in np.load("Na_ops.npy")
        push!(Na_ops, PauliMask(PauliString(i)))
    end
    for i in np.load("Nb_ops.npy")
        push!(Nb_ops, PauliMask(PauliString(i)))
    end
else
    for i in np.load(dir*"ham_ops.npy")
        push!(ham_ops, PauliString(i))
    end
    for i in np.load(dir*"ansatz_ops.npy")
        push!(ansatz_ops, PauliString(i))
    end
    for i in np.load("Na_ops.npy") 
        push!(Na_ops, PauliString(i))
    end
    for i in np.load("Nb_ops.npy") 
        push!(Nb_ops, PauliString(i))
    end
end
