#under construction
using DifferentialEquations
using CairoMakie
using StaticArrays
using SphericalHarmonicExpansions

##
@polyvar x y z # define polynomial variables
# If you have your own fitted polynomical coefficients from Andris's codes, replace the following two lines with your own fitted polynomials by running:
# RF_c=SphericalHarmonicCoefficients(C)
# RF_p = sphericalHarmonicsExpansion(c,x,y,z)
RF_p=0.1*(x^2 - y^2)  # Replace it with your own fitted polynomial
DC_p=-0.001*0.5*(x^2 + y^2 - 2z^2)    # Replace it with your own fitted polynomial

RF_fx,RF_fy,RF_fz=-differentiate(RF_p,(x,y,z))
DC_fx,DC_fy,DC_fz=-differentiate(DC_p,(x,y,z))
myfastfunc(p)=p==0 ? x->0.0 : fastfunc(p)

RF_fx_func=myfastfunc(RF_fx)
RF_fy_func=myfastfunc(RF_fy)
RF_fz_func=myfastfunc(RF_fz)
DC_fx_func=myfastfunc(DC_fx)
DC_fy_func=myfastfunc(DC_fy)
DC_fz_func=myfastfunc(DC_fz)


##
function oscillator(u, p, t)
    # t is the time variable, p is the parameters in named tuple form
    vx,vy,vz,x,y,z= u # unpacking the state vector to  velocity and position
    # v = @view u[1:3]
    r = @view u[4:6]
    # make sure RF_f_func and DC_f_func now return values with respect to the right units
    time_dependent_factor=cos(p.Ω*t)
    dvx=p.RF_fx_func(r)*time_dependent_factor+p.DC_fx_func(r)
    dvy=p.RF_fy_func(r)*time_dependent_factor+p.DC_fy_func(r)
    dvz=p.RF_fz_func(r)*time_dependent_factor+p.DC_fz_func(r)
    dx=vx
    dy=vy
    dz=vz
    SA[dvx,dvy,dvz,dx,dy,dz] # out-of-place version using StaticArrays to make if even faster than in-place version when the number of equations is small
end


##
# define physical parameters
tspan=(0.0,200.0) # time span of the simulation, from 0.0 to 10.0 time units
saveat=0.1 # save the solution every 0.1 time units
u0=SA[0.0,0.0,0.0,1.0,1.0,1.0] # initial conditions, velocity is 0.0, position is 0.1
p=(RF_fx_func=RF_fx_func,RF_fy_func=RF_fy_func,RF_fz_func=RF_fz_func,
   DC_fx_func=DC_fx_func,DC_fy_func=DC_fy_func,DC_fz_func=DC_fz_func,
   Ω=2) # parameters (functions) generated from the fitted polynomials
prob=ODEProblem(oscillator,u0,tspan,p) # define the problem


##
# define numerical parameters
alg=Tsit5() # Tsit5 is Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).
#alg=Vern9() # Vern9 is Verner's “Most Efficient” 9/8 Runge-Kutta method. (lazy 9th order interpolant). It is more efficent when higher accuracy is needed.
abstol=1e-8 # absolute tolerance
reltol=1e-6 # relative tolerance
# smaller tolerances will usually result in more accurate solutions, but will take longer to compute
maxiters=1e5 # maximum number of iterations, this is useful to prevent infinite loops in case of a bug, but you need to set it high enough for the problem to be solved
progress=true # show progress bar or not


##
sol=solve(prob,alg,abstol=abstol,reltol=reltol,saveat=saveat,progress=progress) # solve the problem


##
fig=Figure(size=(900,400))
ax=Axis(fig[1,1],title="Mathieu Oscillator",xlabel="Time",ylabel="x velocity")
lines!(ax,sol.t,sol[1,:],label="x velocity",color=:blue)
ax2=Axis(fig[1,2],title="Mathieu Oscillator",xlabel="Time",ylabel="x position")
lines!(ax2,sol.t,sol[4,:],label="x position",color=:red)
fig


##
# Benchmark to see the allocations and time taken to solve the problem
# If the function is OK, the allocations should be low and not scale with tspan.
using BenchmarkTools
@benchmark solve($prob,$alg,abstol=$abstol,reltol=$reltol,save_on=false) # solve the problem