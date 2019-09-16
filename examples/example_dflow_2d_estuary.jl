# Test case with very simplified 2d estuary (essentially a channel with tidal boundary)
# Reads flow data from netcdf map ouput files from delft3d-fm

using Particles
#include("particles.jl") 
#include("dflow.jl")

#collected configuration is in Dict d 
d=default_userdata() # start with some defaults
#settings for this experiment
n=30 #number of particles
d["nparticles"]=n
# all variables for one particle are collected in a vector 
variables=["x","y","z","age"]
d["variables"]=variables
# initial position of the particles
m=length(variables)
p=zeros(m,n)
p[1,:]=25000.0.+50.0*randn(n,1)
p[2,:]=250.0.+50.0*randn(n,1)
d["particles"]=p #initial values
# simulation time
d["dt"]=3600.0
d["tstart"]=0.0
d["tend"]=25*3600.0
d["plot_map_times"]=collect(0.0:3600.0:(25*3600.0))
d["plot_maps"]=true

#get flow interpolation functions
datadir="test_data"
if !isdir("test_data")
   datadir="../test_data"
end
dflow_map=load_nc_info(datadir,r"estuary_...._map.nc")
interp=load_dflow_grid(dflow_map,50,false)
t0=get_reftime(dflow_map)
u,v=initialize_interpolation(dflow_map,interp,t0)


"""
   !f(ds,s,t,i,d)

Dynamic model, computes time derivative ds of s at current time t
for particle i and possibly using data/functions from d of type userdata.
"""
function f!(dt,s,t,i,d)
   x,y,z,age = s
   # dx/dt=u
   dt[1] = dx   = u(x,y,z,t)
   # dy/dt=v
   dt[2] = dy   = v(x,y,z,t)
   # dz/dt=0
   dt[3] = dz   = 0.0
   # age=(t-t0)
   dt[4] = dage = 1.0
end
d["f"]=f!

function plot_background(d)
   x1=[0.0,100000.0,100000.0,0.0,0.0]
   y1=[0.0,0.0,500.0,500.0,0.0]
   f=plot(x1,y1,legend=false)
   return(f)
end
d["plot_maps_background"]=plot_background

nothing
