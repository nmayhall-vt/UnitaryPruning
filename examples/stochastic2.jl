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

    omat1 = U1'*omat*U1
    omat2 = cos(ga) .* omat .- 1im * sin(ga) .* gmat*omat
    
    println()
    println("transformation expression:")
    display(norm(omat1 - omat2))

    ket = zeros(Bool,N)
    e1 = expectation_value_sign(o,ket)
    e2 = expectation_value_sign(multiply(g,o),ket)

    println()
    println("numerical")
    display(omat1[1])
    println()
    println("analytical")
    display(cos(ga)*e1 - 1im * sin(ga)*e2)
    
end

run()


# vals = [];
# for i in 1:40
#     push!(vals, run(α=0.05, k=i))
# end

# plot(vals, marker = :circle)