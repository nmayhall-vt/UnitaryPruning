using Base
using UnitaryPruning
using Test
using BenchmarkTools 
using PauliOperators
using LinearAlgebra
using Printf
using DifferentialEquations
using OrdinaryDiffEq
using Plots 
using SparseArrays

# @testset "heisenberg" begin
function PauliOperators.expectation_value(ps::PauliSum{N}, ket::KetBitString{N}) where N
    expval = 0
    for (oi,coeff) in ps.ops
        expval += coeff*expectation_value(oi, ket)
    end
    return expval
end

function Base.one(ps::PauliSum{N}) where N 
    return PauliSum(Pauli(0,0,N))
end
function Base.oneunit(ps::PauliSum{N}) where N 
    return PauliSum(Pauli(0,0,N))
end

function SparseArrays.sparse(ps::PauliSum{N}) where N
    a = spzeros(ComplexF64, 2^N, 2^N)
    for (op,coeff) in ps.ops
        a[op.pauli.z,op.pauli.x] = get_coeff
    end
    return a
end


"""
    commutator(a::PauliSum{N}, b::PauliSum{N}) where N

TBW
"""
function commutator(a::PauliSum{N}, b::PauliSum{N}) where N
    out = PauliSum(N)
    for (op_a, coeff_a) in a.ops
        for (op_b, coeff_b) in b.ops
            commute(op_a, op_b) == false || continue
            out += 2*(coeff_a * op_a) * (coeff_b * op_b) 
        end
    end
    return out
end

function heis()
    N = 20 
    ket = KetBitString(N, 0)
    gi, ci = UnitaryPruning.heisenberg(N, 1, 1, 0)
    o = Pauli(N, Z=[1])
    ot = PauliSum(o)

    stepsize = .1
    for ti in 1:100
        ot = UnitaryPruning.bfs_evolution(gi, ci*stepsize, ot, thresh=1e-5)
        l2 = 0
        # for i in Vd
        po = 0
        if haskey(ot,o.pauli) po = ot[o] end
        ev = expectation_value(ot,ket)
        @printf(" t = %12.8f pO = %12.8f ev = %12.8f N = %i\n", ti*stepsize, po, ev, length(ot))
    end
    @show ci
    @show gi
end
function heis_pauli()
    N = 4 
    ket = KetBitString(N, 0)
    gi, ci = UnitaryPruning.heisenberg(N, 1, .1, 0)
    o = Pauli(N, Z=[1])
    ot = PauliSum(o)
    
    display(oneunit(ot))
    display(typeof(one(ot)))
    display(typeof(oneunit(ot)))
    # return
    
    # Now solve DiffEq
    H = PauliSum(N)
    for i in 1:length(gi)
        H += ci[i]*gi[i]
    end
    if is_hermitian(H) == false
        throw(error())
    end

    function time_evolution(ot, p, t)
        # odot = 1im * commutator(H,ot)
        # clip!(odot, thresh=1e-3)
        return ot
    end
    
    tspan = (0.0,10.0)
    prob = ODEProblem(time_evolution,ot,tspan)

    # sol = solve(prob,Tsit5())
    sol = solve(prob)
    # v =[sol.u[i][1] for i in 1:length(sol.u)] 
    # # v =[vmat'*sol.u[i]*vmat for i in 1:length(sol.u)]
    
    # plot(sol.t, real.(v), show=true)
    # print(real.(v))
    # # plot(v,sol.t)

    # U = exp(-1im*10*Hmat)

    # println()
    # display(vmat'*U'*ot*U*vmat)
end

function heis_mat()
    N = 4 
    ket = KetBitString(N, 0)
    gi, ci = UnitaryPruning.heisenberg(N, 1, .1, 0)
    o = Pauli(N, Z=[1])
    ot = PauliSum(o)
    
    # Now solve DiffEq
    H = PauliSum(N)
    for i in 1:length(gi)
        H += ci[i]*gi[i]
    end
    if is_hermitian(H) == false
        throw(error())
    end

    Hmat = Matrix(H)
    vmat = Vector(ket)
    ot = Matrix(ot)
    
    function time_evolution(ot, p, t)
        odot = 1im * (Hmat*ot - ot*Hmat)
        # clip!(odot, thresh=1e-3)
        return odot
    end

    # ot = PauliSum(o)
    # function time_evolution(odot,ot,t)
    #     odot = 1im * commutator(H,ot)
    #     clip!(odot, thresh=1e-3)
    #     return odot
    # end
    tspan = (0.0,10.0)
    prob = ODEProblem(time_evolution,ot,tspan)

    # sol = solve(prob,Tsit5())
    sol = solve(prob)
    v =[sol.u[i][1] for i in 1:length(sol.u)] 
    # v =[vmat'*sol.u[i]*vmat for i in 1:length(sol.u)]
    
    plot(sol.t, real.(v), show=true)
    print(real.(v))
    # plot(v,sol.t)

    U = exp(-1im*10*Hmat)

    println()
    display(vmat'*U'*ot*U*vmat)
end

# heis()
heis_pauli()
# heis_mat()