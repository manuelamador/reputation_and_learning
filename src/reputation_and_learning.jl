using Roots
using Parameters
using DifferentialEquations
using QuadGK


const _tmax = 1000.0

abstract type AbstractModel end 


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
    @unpack r, i, λ, bmax = m
    return max((r  - (i + λ * (1  - q)) / q), 0) * (bmax - b)
end 

function bigQ(m::SimpleModel, b, c)
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
    @unpack i, y, λ, r = m
    (r * (c + i - y + λ))/(r - b * r + λ)^2
end 



function cF(m::AbstractModel, b, q)
    @unpack r, i, λ, y = m

    return y - (i + λ) * b + q * (h(m, b, q) + λ * b)
end 


function bPrime(m::AbstractModel, b, c)
    h(m, b, bigQ(m, b, c))
end 


function qPrime(m::AbstractModel, b, c)
    bigQ_prime(m, b, c) * h(m, b, bigQ(m, b, c))
end 


function x(m::AbstractModel, b, c)
    @unpack i, λ = m 
    (qPrime(m, b, c) + (i + λ) * (1 - bigQ(m, b, c))) / bigQ(m, b, c)
end 


function ρ_prime(m::AbstractModel, ρ, b, c)
    @unpack ϵ, δ = m
    (1 - ρ) * ϵ + ρ * (x(m, b, c) - δ)
end 


function solve_b_ODE(m::AbstractModel, c; tmax=_tmax)
    condition(b, t, integrator) = bigQ(m, b, c) - qupperbar(m)
    affect!(integrator) = terminate!(integrator)

    cb = ContinuousCallback(condition, affect!)

    tspan = (0.0, tmax)
    prob = ODEProblem((u, p, t) -> bPrime(m, u, c), 0.0, tspan)

    solve(prob, Tsit5(), callback=cb)
end 


function solve_ρ_ODE(m::AbstractModel, c, b_sol)
    tspan = (0.0, b_sol.t[end])
    prob = ODEProblem(
        (ρ, p, t) -> ρ_prime(m, ρ, b_sol(t), c), 0.0, tspan)
    solve(prob, Tsit5())
end 


function check_sol(m::AbstractModel, c)
    bsol = solve_b_ODE(m, c)
    ρsol = solve_ρ_ODE(m, c, bsol)
    ρsol.u[end] - 1
end 


function find_c_star(m::AbstractModel)
    find_zero(c -> check_sol(m,c), (1.0000001, cbar(m)))
end 


function solve_b_after_T(m::AbstractModel, b_sol, T; tmax=_tmax)
    tspan = (T, _tmax)
    prob = ODEProblem(
        (b, p, t) -> h(m, b, qupperbar(m)), b_sol.u[end], tspan)
    solve(prob, Tsit5())
end


function solve_F_ODE(m::AbstractModel, c, bsol, ρsol)
    tspan = (0.0, bsol.t[end])
    prob = ODEProblem(
        (F, p, t) -> x(m, bsol(t), c) * (1 - F) / (1 - ρsol(t)), 0.0, tspan)
    solve(prob, Tsit5())
end 


function construct_solution(m::AbstractModel; tmax=_tmax)
    @unpack i, λ, y, δ = m

    cstar = find_c_star(m)
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
    trade_deficit=trade_deficit, x=xF)
end


function plots(sol; plot_tmax=80)
    trange = collect(range(0, plot_tmax, length=100))

    p1 = plot(trange, [sol.q(t) for t ∈ trange], legend=false, lw=2, 
        title="q")
    vline!(p1, [sol.T])


    p2 = plot(trange, [sol.b(t) for t ∈ trange], legend=false, lw=2, 
    title="b")
    vline!(p2, [sol.T])

    p3 = plot(trange, [sol.ρ(t) for t ∈ trange], legend=false, lw=2, 
    title="ρ")
    vline!(p3, [sol.T])

    p4 = plot(trange, [sol.c(t) for t ∈ trange], legend=false, lw=2, 
    title="c")
    vline!(p4, [sol.T])

    p5 = plot(trange, [sol.yield(t) for t ∈ trange], legend=false, lw=2, 
    title="yield")
    vline!(p5, [sol.T])
    
    p6 = plot(trange, [sol.F(t) for t ∈ trange], legend=false, lw=2, 
    title="\$F_0\$")
    vline!(p6, [sol.T])
    
    p7 = plot(trange, [sol.default_rate(t) for t ∈ trange], legend=false, lw=2, 
    title="default rate")
    vline!(p7, [sol.T])

    p8 = plot(trange, [sol.trade_deficit(t) for t ∈ trange], legend=false, lw=2, 
    title="trade deficit")
    vline!(p8, [sol.T])

    p9 = plot(trange, [sol.x(t) for t ∈ trange], legend=false, lw=2, 
    title="x")
    vline!(p9, [sol.T])

    plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=(3, 3), size=(1000, 800))
end 


function plots(m::AbstractModel; tmax=_tmax, plot_tmax=80)

    sol = construct_solution(m, tmax=tmax)
    plots(sol, plot_tmax=plot_tmax)
end 


function get_m(sol) 
    T = sol.T
    x = sol.x 
    f = t -> exp(quadgk(x, t, T)[1])
    return T + quadgk(t -> (t * x(t) * f(t)), 0, T)[1]
end 