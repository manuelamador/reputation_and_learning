#  Basic type and methods to solve for the results in the "Reputation and 
#  Sovereign Default" paper by Manuel Amador and Christopher Phelan. 
# 

using Roots
using Parameters
using DifferentialEquations
using QuadGK
using Plots
using Interpolations

# The maximum range of t for the solutions of the ODE. 
const _tmax = 100.0

abstract type AbstractModel end 

# The particular simplemodel that uses the h function described in the paper. 
@with_kw struct SimpleModel{T <: Real} <: AbstractModel
    r::T = 0.15
    i::T = 0.01
    λ::T = 0.2
    δ::T = 0.02
    ϵ::T = 0.01
    y::T = 1.0
    bmax::T = y
    qupperbar::T = (i + λ)/(i + λ + δ) 
    qlowerbar::T = (i + λ)/(r + λ)
    cbar::T = r - i + y
end 

function h(m::SimpleModel, b, q)
    # The function h assumed in the paper. 
    @unpack r, i, λ, bmax = m
    return max((r  - (i + λ * (1  - q)) / q), 0) * (bmax - b)
end 

function bigQ(m::SimpleModel, b, c)
    # bigQ is obtained by solving the equation cF(b, q) = c. 
    @unpack i, y, λ, r = m
    (c + i - y + λ)/(r - b * r + λ)
end 

function qlowerbar(m::SimpleModel)
    m.qlowerbar
end 

function qupperbar(m::SimpleModel)
    m.qupperbar 
end 

function cbar(m::SimpleModel)
    m.cbar
end 


function bigQ_prime(m::SimpleModel, b, c)
    # The analytical derivative of Q.
    @unpack i, y, λ, r = m
    (r * (c + i - y + λ))/(r - b * r + λ)^2
end 


function cF(m::AbstractModel, b, q)
    # The consumption function. 
    @unpack r, i, λ, y = m
    return y - (i + λ) * b + q * (h(m, b, q) + λ * b)
end 


function bPrime(m::AbstractModel, b, c)
    # The time derivative of b given a constant consumption level c. 
    h(m, b, bigQ(m, b, c))
end 


function qPrime(m::AbstractModel, b, c)
    # The (time) derivative of Q given constant consumption c. 
    bigQ_prime(m, b, c) * h(m, b, bigQ(m, b, c))
end 


function x(m::AbstractModel, b, c)
    # The implied default rate given a constant consumption c. 
    @unpack i, λ = m 
    (qPrime(m, b, c) + (i + λ) * (1 - bigQ(m, b, c))) / bigQ(m, b, c)
end 


function ρ_prime(m::AbstractModel, ρ, b, c)
    # The time derivative of market belief ρ given constant consumption c. 
    @unpack ϵ, δ = m
    (1 - ρ) * ϵ + ρ * (x(m, b, c) - δ)
end 


function solve_b_ODE(m::AbstractModel, c; tmax=_tmax)
    # Given an initial consumption, c, solves the evolution of debt, b(t)
    # up to the point where the price Q riches qupperbar
    condition(b, t, integrator) = bigQ(m, b, c) - qupperbar(m)
    affect!(integrator) = terminate!(integrator)

    cb = ContinuousCallback(condition, affect!)

    tspan = (0.0, tmax)
    prob = ODEProblem((u, p, t) -> bPrime(m, u, c), 0.0, tspan)

    solve(prob, Tsit5(), callback=cb)
end 


function solve_ρ_ODE(m::AbstractModel, c, b_sol)
    # Given a solution for b and an initial condition, compute the 
    # evolution for ρ.
    tspan = (0.0, b_sol.t[end])
    prob = ODEProblem(
        (ρ, p, t) -> ρ_prime(m, ρ, b_sol(t), c), 0.0, tspan)
    solve(prob, Tsit5())
end 


function eqm_distance(m::AbstractModel, c0)
    # Solve the model given initial c0 and returns the difference 
    # between ρ and 1 at cap T. In eqm, ρ = 1. 
    bsol = solve_b_ODE(m, c0)
    ρsol = solve_ρ_ODE(m, c0, bsol)
    ρsol.u[end] - 1
end 


function find_cstar(m::AbstractModel)
    find_zero(c -> eqm_distance(m,c), (1.0000001, cbar(m)))
end 


function solve_b_after_T(m::AbstractModel, b_sol, T; tmax=_tmax)
    # The evolution of debt after T. The price is at qupperbar. 
    tspan = (T, _tmax)
    prob = ODEProblem(
        (b, p, t) -> h(m, b, qupperbar(m)), b_sol.u[end], tspan)
    solve(prob, Tsit5())
end


function solve_F_ODE(m::AbstractModel, c, bsol, ρsol)
    # Integrating to obtain the F (default CDF) function. 
    tspan = (0.0, bsol.t[end])
    prob = ODEProblem(
        (F, p, t) -> x(m, bsol(t), c) * (1 - F) / (1 - ρsol(t)), 0.0, tspan)
    solve(prob, Tsit5())
end 


function construct_solution(m::AbstractModel; tmax=_tmax)
    # Find the initial c0, and construct the equilibrium objects used to 
    # make the figures. 
    @unpack i, λ, y, δ = m

    cstar = find_cstar(m)
    bsol = solve_b_ODE(m, cstar, tmax=tmax)
    ρsol = solve_ρ_ODE(m, cstar, bsol)    
    Fsol = solve_F_ODE(m, cstar, bsol, ρsol)

    T = bsol.t[end]
    bsol_after_T = solve_b_after_T(m, bsol, T, tmax=tmax)

    
    q = t -> (t <= T ? bigQ(m, bsol(t), cstar) : qupperbar(m))
    b = t -> (t <= T ? bsol(t) : bsol_after_T(t))
    ρ = t -> (t <= T ? ρsol(t) : 1.0)
    c = t -> cF(m, b(t), q(t))
    yield = t ->  (i + λ) / q(t) - λ
    F = t -> (t <= T ? Fsol(t) : 1.0) 
    default_rate = t -> (t <= T ? x(m, b(t), cstar) / (1 - ρ(t)) : Inf)
    trade_deficit = t -> 100 * (c(t) - y) / y 
    xF = t -> (t <= T ? x(m, b(t), cstar) : δ)


    (bsol=bsol, ρsol=ρsol, T=T, q=q, b=b, ρ=ρ, c=c, yield=yield, F=F, default_rate=default_rate, 
    trade_deficit=trade_deficit, x=xF, m=m)
end


function figure1(sol; plot_tmax=80)
    trange = collect(range(0, plot_tmax, length=500))
    plist = []
    for (var, name) in ( 
                        (sol.q, "q"), 
                        (sol.b, "b"),
                        (sol.ρ, "ρ"), 
                        (sol.c, "c"),
                        (sol.yield, "yield"),
                        (sol.F, "\$F_0\$"),
                        (sol.default_rate, "default rate"),
                        (sol.trade_deficit, "trade deficit"),
                        (sol.x, "x"))          
        p = plot(trange, [var(t) for t ∈ trange], 
                legend=false, lw=4, 
                title=name)
        vline!(p, [sol.T])
        push!(plist, p)
    end
    
    plot(plist..., layout=(3, 3), size=(1000, 800))
end 


function figure1(m::AbstractModel; tmax=_tmax, plot_tmax=80)
    # If we dont have the solution, compute first
    sol = construct_solution(m, tmax=tmax)
    # then do the figure. 
    figure1(sol, plot_tmax=plot_tmax)
end 


function get_m(sol) 
    # Computes the m value by solving the double integral stated in the paper. 
    T, x = sol.T, sol.x 
    f = t -> exp(quadgk(x, t, T)[1])
    return T + quadgk(t -> (t * x(t) * f(t)), 0, T)[1]
end 


function figure2(sol)
    @unpack T  = sol
    bF, qF, ρF = sol.b, sol.q, sol.ρ 
    @unpack ϵ, δ, i, λ = sol.m

    trange = range(0.0, T, length=100)

    # Creating ρ as a function of b
    ρ_interp =LinearInterpolation( [bF(t) for t in trange], [ρF(t) for t in trange])
    # Creating q as a function of b
    q_interp =LinearInterpolation( [bF(t) for t in trange], [qF(t) for t in trange])

    brange = range(0.0, bF(T), length=17)[2:end-1]

    # We are going to solve backwards a differential system with two states ρ and q. 
    # where both ρ and q are functions of b. 
    #
    # That means that, for the differential equations package, 
    # u = [ρ, q]
    # t  = b    

    condition(u, t, integrator) = u[1] - 1 # Stop integration if ρ reaches 1. 
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    sols = []
    for bi in brange 
        ρ0 = ρ_interp(bi)  # initial condition for ρ
        q0 = q_interp(bi)  # initial condition for q

        # The range of b (recall, we are solving backwards)
        bspan = (bi, 0.0)

        # Recall: 
        # u = [ρ, q]
        # t  = b 
        #   
        # The ODE system for ρ and q as functions of b is given by: 
        ρb_prime = (ρ, q, b) ->  (ϵ * (1 - ρ) - ρ * δ) / h(m, b, q)
        qb_prime = (q, b) -> ((i + λ) * q - (i + λ)) / h(m, b, q)
        prob = ODEProblem(
            (u, p, t) -> [ρb_prime(u[1], u[2], t), qb_prime(u[2], t)],
            [ρ0, q0],
            bspan
        )
        push!(sols, solve(prob, Tsit5(), callback=cb)) # store the solution.
    end 

    p = plot()
    for s in sols 
        plot!(p, s.t, [s(b)[1] for b in s.t], lw=2, color=:gray)
    end 
    plot!(p, [bF(t) for t in trange], [ρF(t) for t in trange], 
            lw=4, 
            color=:blue, 
            legend=false, 
            size=(600, 400))
    p
end


function figure2(m::AbstractModel; tmax=_tmax)
    # If we don't have the solution, compute it first:
    sol = construct_solution(m, tmax=tmax)
    # then do the plot
    figure2(sol)
end 