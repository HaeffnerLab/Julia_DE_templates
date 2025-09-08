#07/10/2025 by Baiyi
using DifferentialEquations
using CairoMakie
using StaticArrays


##
function harmonic_oscillator(u, p, t)
    # t is the time variable, p is the parameters in named tuple form
    v,x= u # unpacking the state vector to v and x, velocity and position
    # dv=-p.ωsq*x*sin(t)- (2pi*1.0)^2/5*x^2# p. ωsq is ω^2, precomputed for efficiency 
    dv=(-p.ωsq*x- (2pi*1.0)^2/5*x^2)*sin(p.ωsq*10*t)
    dx=v
    SA[dv,dx] # out-of-place version using StaticArrays to make if even faster than in-place version when the number of equations is small
end


##
# define physical parameters
ω=2pi*1.0
tspan=(0.0,1000.0) # time span of the simulation, from 0.0 to 10.0 time units
saveat=0.1 # save the solution every 0.1 time units
u0=SA[0.0,20.0] # initial conditions, velocity is 0.0, position is 1.0
p=(ωsq=ω^2,) # parameters, ω^2 is precomputed for efficiency
prob=ODEProblem(harmonic_oscillator,u0,tspan,p) # define the problem

# Create x range for energy plot
x_range = range(-30, 30, length=1000)
# Calculate energy values
E_values = [0.5*p.ωsq*x^2+ (2pi*1.0)^2/5*x^3/3 for x in x_range]

# Plot the energy function
fig_energy = Figure(size=(600,400))
ax_energy = Axis(fig_energy[1,1], title="Energy Function", xlabel="Position (x)", ylabel="Energy (E)")
lines!(ax_energy, x_range, E_values, label="E = 0.5ω²x² + (2π)²/5 * x³/3", color=:red, linewidth=2)

# Add a point at the initial position x = u0[2] (which is -7.0)
initial_x = u0[2]  # This is -7.0 from your initial conditions
initial_E = 0.5*p.ωsq*initial_x^2 + (2pi*1.0)^2/5*initial_x^3/3
scatter!(ax_energy, [initial_x], [initial_E], color=:blue, markersize=10, label="Initial position")

axislegend(ax_energy)
fig_energy


##
# define numerical parameters
alg=Vern9() # Tsit5 is Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).
#alg=Vern9() # Vern9 is Verner's “Most Efficient” 9/8 Runge-Kutta method. (lazy 9th order interpolant). It is more efficent when higher accuracy is needed.
abstol=1e-6 # absolute tolerance
reltol=1e-4 # relative tolerance
# smaller tolerances will usually result in more accurate solutions, but will take longer to compute
maxiters=1e5 # maximum number of iterations, this is useful to prevent infinite loops in case of a bug, but you need to set it high enough for the problem to be solved
progress=true # show progress bar or not


##
sol=solve(prob,alg,abstol=abstol,reltol=reltol,saveat=saveat,progress=progress) # solve the problem
println("Solution computed successfully!")


fig=Figure(size=(900,400))
ax=Axis(fig[1,1],title="Harmonic Oscillator",xlabel="Time",ylabel="velocity")
lines!(ax,sol.t,sol[1,:],label="velocity",color=:blue)
ax2=Axis(fig[1,2],title="Harmonic Oscillator",xlabel="Time",ylabel="position")
lines!(ax2,sol.t,sol[2,:],label="position",color=:red)
println("Figure created successfully!")
fig

##
# Benchmark to see the allocations and time taken to solve the problem
# If the function is OK, the allocations should be low and not scale with tspan.
using BenchmarkTools
@benchmark solve($prob,$alg,abstol=$abstol,reltol=$reltol,save_on=false) # solve the problem
# println("Benchmark completed successfully!")