# Test case with very simplified 2d estuary (essentially a channel with tidal boundary)
# Reads flow data from netcdf map ouput files from delft3d-fm

using Particles
#include("particles.jl") 
#include("dflow.jl")
#include("wms_client.jl")

#collected configuration is in Dict d 
d=default_userdata() # start with some defaults
#settings for this experiment
n=30 #number of particles
d["nparticles"]=n
# all variables for one particle are collected in a vector 
d["coordinates"]="spherical" #projected or spherical 
variables=["lon","lat","age"]
d["variables"]=variables
d["bbox"]=[-15.0,13.0,64.0]
# initial position of the particles
m=length(variables)
p=zeros(m,n)
p[1,:]=5.0.+0.05*randn(n,1)
p[2,:]=53.5.+0.05*randn(n,1)
d["particles"]=p #initial values
# simulation time
# reference 22-12-2012 + 240-408 hours
d["reftime"]=DateTime(2012,12,22) 
h=3600.0 #seconds per hour
d["dt"]=600.0 #seconds
d["tstart"]=240.0*h
d["tend"]=408.0*h
#write to netcdf
d["write_maps_times"]=collect((240.0*h):(1.0*h):(408*h))
d["write_maps"]=true
d["write_maps_filename"]="output_dflow_2d_dcsm.nc"
#plot to screen 
d["plot_map_times"]=collect((240.0*h):(1.0*h):(408*h))
d["plot_maps"]=false

#get flow interpolation functions
datadir="test_data"
if !isdir(datadir)
   datadir="../test_data"
end
dflow_map=load_nc_info(datadir,r"DCSM-FM_0_5nm_...._map.nc")
interp=load_dflow_grid(dflow_map,50,false)
t0=get_reftime(dflow_map)
u,v=initialize_interpolation(dflow_map,interp,t0)


"""
   !f(ds,s,t,i,d)

Dynamic model, computes time derivative ds of s at current time t
for particle i and possibly using data/functions from d of type userdata.
"""
function f!(dt,s,t,i,d)
   lon,lat,age = s
   z=0.0
   R=6371.0e3 #mean radius of earth from wikipedia
   deg2rad=pi/180.0 #convert degrees to radians
   # dx/dt=u
   dt[1] = dlon = u(lon,lat,z,t)/(R*cos(deg2rad*lat))
   # dy/dt=v
   dt[2] = dlat = v(lon,lat,z,t)/R
   # age=(t-t0)
   dt[3] = dage = 1.0
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
