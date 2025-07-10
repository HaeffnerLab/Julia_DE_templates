using DifferentialEquations
using CairoMakie
using StaticArrays


##
function harmonic_oscillator(u, p, t)
    # t is the time variable, p is the parameters in named tuple form
    v,x= u # unpacking the state vector to v and x, velocity and position
    dv=-p.ωsq*x # p. ωsq is ω^2, precomputed for efficiency 
    dx=v
    SA[dv,dx] # out-of-place version using StaticArrays to make if even faster than in-place version when the number of equations is small
end


##
# define physical parameters
ω=2pi*1.0
tspan=(0.0,10.0) # time span of the simulation, from 0.0 to 10.0 time units
saveat=0.1 # save the solution every 0.1 time units
u0=SA[0.0,1.0] # initial conditions, velocity is 0.0, position is 1.0
p=(ωsq=ω^2,) # parameters, ω^2 is precomputed for efficiency
prob=ODEProblem(harmonic_oscillator,u0,tspan,p) # define the problem


##
# define numerical parameters
alg=Tsit5() # Tsit5 is Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).
#alg=Vern9() # Vern9 is Verner's “Most Efficient” 9/8 Runge-Kutta method. (lazy 9th order interpolant). It is more efficent when higher accuracy is needed.
abstol=1e-6 # absolute tolerance
reltol=1e-4 # relative tolerance
# smaller tolerances will usually result in more accurate solutions, but will take longer to compute
maxiters=1e5 # maximum number of iterations, this is useful to prevent infinite loops in case of a bug, but you need to set it high enough for the problem to be solved
progress=true # show progress bar or not


##
sol=solve(prob,alg,abstol=abstol,reltol=reltol,saveat=saveat,progress=progress) # solve the problem


##
fig=Figure(size=(900,400))
ax=Axis(fig[1,1],title="Harmonic Oscillator",xlabel="Time",ylabel="velocity")
lines!(ax,sol.t,sol[1,:],label="velocity",color=:blue)
ax2=Axis(fig[1,2],title="Harmonic Oscillator",xlabel="Time",ylabel="position")
lines!(ax2,sol.t,sol[2,:],label="position",color=:red)
fig


##
# Benchmark to see the allocations and time taken to solve the problem
# If the function is OK, the allocations should be low and not scale with tspan.
using BenchmarkTools
@benchmark solve($prob,$alg,abstol=$abstol,reltol=$reltol,save_on=false) # solve the problem