using UnitaryPruning
using Plots
using Statistics
using Printf
using Random
using LinearAlgebra

function run(;ga=1.23)
    N = 8 

    # I = PauliBoolVec(N)
    g = PauliBoolVec(N,X=[1,2], Y=[3], Z=[5,6])
    o = PauliBoolVec(N,X=[3], Y=[1,2], Z=[5,8])
    ga = 1.23

    omat = to_matrix(o)
    gmat = to_matrix(g)
    Imat = to_matrix(PauliBoolVec(N))
    display(commute(o,g))
    
    println()
    println("multiplication:")
    display(norm(omat*gmat - to_matrix(multiply(o,g))))
    display(norm(omat*gmat - to_matrix(multiply(o,g))))
    display(norm(gmat*gmat - Imat))
    display(norm(to_matrix(multiply(g,g)) - Imat))

    # Test exponential
    U1 = exp(1im * ga * gmat ./ 2)
    U2 = (cos(ga/2) .* Imat) .+ (1im * sin(ga/2) .* gmat)
    
    println()
    println("exponential expression:")
    display(norm(abs.(U1)-abs.(U2)))
    display(norm(U1-U2))

    println()
    println("exponential determinants:")
    display(det(U1))
    display(det(U2))
 
    # Test transformation

    omat1 = U1*omat*U1'
    omat2 = cos(ga) .* omat .+ 1im * sin(ga) .* gmat*omat
    
    println()
    println("transformation expression:")
    display(norm(omat1 - omat2))

    ket = zeros(Bool,N)
    e1 = expectation_value_sign(o,ket)
    e2 = expectation_value_sign(multiply(g,o),ket)

    println("numerical")
    display(omat1[1])
    println("analytical")
    display(cos(ga)*e1 + 1im * sin(ga)*e2)
    return
    Nt = length(generators)
    length(angles) == Nt || throw(DimensionMismatch)
    for t in 1:Nt
        α = angles[t]
        # Ut = e(i α Pn) = cos(α) I + i sin(α) Pn
        u = cos(α) .* u .+ 1im*sin(α) .* u * to_matrix(generators[t])
    end
    #Mz
    o = PauliBoolVec(N, Z=[1])
    o_mat = to_matrix(o)
    # for i in 2:N
    #     o .+= to_matrix(PauliBoolVec(N, Z=[i]))
    # end


    U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    
    # Ui = UnitaryPruning.build_time_evolution_matrix([generators[2]], [parameters[2]])
    # Uj = exp(to_matrix(generators[2]) .* 1im .* parameters[2])

    # display(norm(Ui - Uj))

    # out = [mean(results[1:i]) for i in 1:length(results)]
    # plot(real(results))
    m = diag(U*o_mat*U')
    println(" expectation values:")
    display(m[1])
    
    ket = zeros(Bool, N)
    Random.seed!(1)
    results = Vector{ComplexF64}([])
    for i in 1:10000
        e, path = UnitaryPruning.stochastic_pauli_dynamics_run(generators, parameters, o, ket)
        # display(path)
        push!(results, e)
    end
    @printf(" Mean: %12.8f Stdev: %12.8f\n", mean(results), std(results))

    out = [results[1]]
    for (idx,i) in enumerate(results)
        idx > 1 || continue
        push!(out, (out[idx-1]*(idx-1)+i)/idx)
    end

    return real(m), out 
end

run()


# vals = [];
# for i in 1:40
#     push!(vals, run(α=0.05, k=i))
# end

# plot(vals, marker = :circle)