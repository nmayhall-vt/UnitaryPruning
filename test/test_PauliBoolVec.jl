using UnitaryPruning
using Test
using BenchmarkTools 

@testset "PauliBoolVec" begin

    a = PauliBoolVec("XYZIXY")
    b = PauliBoolVec("YXXYZZ")
    c = PauliBoolVec("ZZYYYX")
    c.θ = 2
    display(c)
    display(a*b)
    @test c == a*b 

    @test commute(a,b) == false
    
    c.θ = 0
    @test c == b*a 

    a = PauliBoolVec("XYYZIZYZIIXYIIYZI")
    b = PauliBoolVec("YYYIIZIIZIZYXIXYZ")
    c = PauliBoolVec("ZIIZIIYZZIYIXIZXZ")
    c.θ = 0
    display(a)
    display(b)
    display(c)
    display(a*b)
    @test c == a*b 
    @test c == b*a
    
    @test commute(a,b) == true 

    display(UnitaryPruning.is_diagonal(a))

    
    a = PauliBoolVec("IIZIZIIII")
    v = Vector{Bool}([1,1,1,0,0,0,0,0,0])
    @test expectation_value_sign(a,v) == -1
    
    v = Vector{Bool}([1,1,0,1,0,0,0,0,0])
    @test expectation_value_sign(a,v) == 1
    
    v = Vector{Bool}([1,1,1,1,1,0,0,0,0])
    @test expectation_value_sign(a,v) == 1




    #############################################################
    g = PauliBoolVec(8,X=[1,2], Y=[3], Z=[5,6])
    o = PauliBoolVec(8,X=[3], Y=[1,2], Z=[5,8])

    gmat = to_matrix(g)
    omat = to_matrix(o)
    Imat = to_matrix(PauliBoolVec(8))
    println()
    println("multiplication:")
    e1 = norm(omat*gmat - to_matrix(multiply(o,g)))
    e2 = norm(gmat*gmat - Imat)
    e3 = norm(to_matrix(multiply(g,g)) - Imat)

    @test abs(e1) < 1e-12
    @test abs(e2) < 1e-12
    @test abs(e3) < 1e-12

    #############################################################


    Random.seed!(1)
    N = 8
    generators = Vector{PauliBoolVec{N}}([])
    angles = rand(20) 
    
    o = UnitaryPruning.random_PauliBoolVec(N)

    for i in 1:20
        push!(generators, UnitaryPruning.random_PauliBoolVec(N))
    end
    U = UnitaryPruning.build_time_evolution_matrix(generators, angles)

    @test abs(abs(det(U))-1) < 1e-12

    U2 = to_matrix(PauliBoolVec(N))
    for t in 1:length(generators)
        U2 = exp(1im * angles[t] .* to_matrix(generators[t]) ./ 2) * U2
    end

    e = norm(U2 - U)
    @test abs(e) < 1e-12

    #############################################################

    o = PauliBoolVec("YZXZYZZY")
    g = PauliBoolVec("XZXZXIZX")
    
    display(o)
    display(g)
    ket = zeros(Bool, N)
    par = 1.1

    println(" Do they commute?: ", commute(o,g))

    U1 = UnitaryPruning.build_time_evolution_matrix([g], [par])
    e0 = to_matrix(o)[1]
    e1 = (U1'*to_matrix(o)*U1)[1]
   
    e2 = expectation_value_sign(o,ket)

    e3list = Vector{ComplexF64}([])
    for i in 1:1000
        ei, _ = UnitaryPruning.stochastic_pauli_dynamics_run([g], [par], o, ket)
        push!(e3list, ei)
    end
  
    plot(real(e3list))
    e3 = mean(e3list)

    @printf(" e0: %12.8f %12.8fi\n", real(e0), imag(e0))
    @printf(" e1: %12.8f %12.8fi\n", real(e1), imag(e1))
    @printf(" e2: %12.8f %12.8fi\n", real(e2), imag(e2))
    @printf(" e3: %12.8f %12.8fi\n", real(e3), imag(e3))

    # display(sin(par))
    # display(1im*sin(par)*expectation_value_sign(multiply(o,g), ket))
end