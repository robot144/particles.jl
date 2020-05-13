# Test case with simple flow in x-direction with vertical coordinate

using Particles
using Plots
#include("particles.jl")

#collected configuration is in Dict d
d=default_userdata() # start with some defaults
#settings for this experiment
n=30 #number of particles
d["nparticles"]=n
# all variables for one particle are collected in a vector
variables=["x","y","z"]
d["variables"]=variables
# initial position of the particles
m=length(variables)
p=zeros(m,n)
p[index("x",variables),:].=0.0
p[index("z",variables),:]=0.0.-10*rand(n,1)
d["particles"]=p #initial values
# simulation time
d["dt"]=0.1
d["tstart"]=0.0
d["tend"]=90.0
#write to netcdf
d["write_maps_times"]=collect(0.0:1.0:90.0)
d["write_maps"]=true
d["write_maps_filename"]="output_2dv_sideview.nc"
#write plots to file
d["plot_maps_times"]=collect(0.0:10.0:90.0)
d["plot_maps"]=true

#flow in x direction (for now has to be called u)
function u(x,y,z,t)
   # u=-s_y
   1.0
end

"""
   !f(ds,s,t,i,d)

Dynamic model, computes time derivative ds of s at current time t
for particle i and possibly using data/functions from d of type userdata.
"""
function f!(dt,s,t,i,d)
   x,y,z = s
   # dx/dt=u
   dt[1] = dx   = u(x,y,z,t)
   # dy/dt=v
   dt[2] = dy   = 0.0
   # dz/dt=w_up floating
   if (z<0.0)
      dt[3] = dz   = 0.15
   else
      dt[3] = dz   = 0.0
   end
end
d["f"]=f!

function plot_background(d)
   plot([0.0,100.0],[-10.0,-10.0])
end
d["plot_maps_background"]=plot_background
d["plot_maps_func"]=plot_maps_xz

println("Start with run_simulation(d) if it does not start automatically")
nothing
