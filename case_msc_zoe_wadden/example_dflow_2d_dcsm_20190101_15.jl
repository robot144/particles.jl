# A test rund with data (DCSM-FM 2D flowdata) for the MSC ZOE accident.
#

using Particles
using Dates
using Plots
#include("particles.jl") 
#include("dflow.jl")
#include("wms_client.jl")

#collected configuration is in Dict d 
d=default_userdata() # start with some defaults
#settings for this experiment
n=1000 #number of particles
d["nparticles"]=n
# all variables for one particle are collected in a vector 
d["coordinates"]="spherical" #projected or spherical 
variables=["lon","lat","age"]
d["variables"]=variables
#d["bbox"]=[-15.0,43.0,13.0,64.0] #model domain
d["bbox"]=[4.5,53.0,5.5,54.0] #where we expect particles
# initial position of the particles
m=length(variables)
p=zeros(m,n)
p[1,:]=5.0.+0.05*randn(n,1)
p[2,:]=53.5.+0.05*randn(n,1)
d["particles"]=p #initial values
# simulation time
# reference 22-12-2018 in file
# times in input file 2019-01-01 - 2019-01-15
d["reftime"]=DateTime(2019,01,01) #we use 2019-01-01 as refdate for this run 
h=3600.0 #seconds per hour
d["dt"]=1800.0 #seconds
d["tstart"]=0.0*h
#d["tend"]=14.0*24.0*h
d["tend"]=2.0*24.0*h
#write to netcdf
d["write_maps_times"]=collect((0.0*h):(1.0*h):(14*24*h)) 
d["write_maps"]=true
d["write_maps_filename"]="output_dflow_2d_dcsm2019.nc"
#plot maps 
d["plot_maps_times"]=collect((0.0*h):(1.0*h):(14*24*h))
d["plot_maps"]=true

"""
   d["f"]=initialie_model(d)
Create model functions with a local scope.
"""
function initialize_model(d)
   datadir="."
   if !isdir(datadir)
      datadir="../testdata"
   end
   #get flow interpolation functions
   dflow_map=load_nc_info(datadir,r"DCSM-FM_0_5nm_...._map.nc")
   interp=load_dflow_grid(dflow_map,50,false)
   #t0=get_reftime(dflow_map)
   t0=d["reftime"]
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
      rad2deg=180.0/pi
      # dx/dt=u
      dt[1] = dlon = rad2deg*u(lon,lat,z,t)/(R*cos(deg2rad*lat))
      # dy/dt=v
      dt[2] = dlat = rad2deg*v(lon,lat,z,t)/R
      # age=(t-t0)
      dt[3] = dage = 1.0
   end
   return f!
end
d["f"]=initialize_model(d)

#prepage background image
plot_maps_size=d["plot_maps_size"]
width,height=plot_maps_size
plot_bbox=d["bbox"]
wms_server=WmsServer("emodnet-bathymetry")
img=get_map(wms_server,plot_bbox,width,height)
d["background_image"]=img
function plot_background(d)
   plot_bbox=d["bbox"]
   img=d["background_image"]
   f=plot_image(img,plot_bbox)
   #f=plot(x1,y1,legend=false)
   return(f)
end
d["plot_maps_background"]=plot_background


println("run_simulation(d) to start run")
nothing
