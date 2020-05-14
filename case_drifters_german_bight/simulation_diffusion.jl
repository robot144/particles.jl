using Particles
using Plots
include("drifterfunctions.jl")

drifternumber = 5;
water_3d = true

d=default_userdata()
n=1801                                                                              #Drfiter pos, using 2Dwwater, using 3d water, 2d jonswap, 3d jonswap, 2d 1.6% wind, 3d 1.6% wind

d["nparticles"]=n
d["coordinates"]="spherical"
# d["bbox"] = [6.0,53.5,9.1,55.5]
d["bbox"] = [9,53.5,12.5,58]                                                  # Where we expect particles
d["plot_maps_size"] = (1500,1000)


# Use different data directories, depending on the drifternumber and water_3d
# Drifter 1-4 are in the period of march/april of 2017, thus dflow and matroos data of this period should be used
# Drifter 5-7 are in the period of october of 2017, thus then this dflow/matroos data should be used
# water_3d is a boolean indicating if 2D or 3D dflow data should be used (true for 3D, false for 2D)
if drifternumber >= 5
   datadir_2d = "data/october_2d"
   datadir_3d = "data/october_3d"
   matroos_data = MatroosData("data","matroos_october.nc")
else
   datadir_2d = "data/march_2d"
   datadir_3d = "data/march_3d"
   matroos_data = MatroosData("data","matroos_march.nc")
end

dflow_map_2d = load_nc_info(datadir_2d, r"DCSM-FM_0_5nm_...._map.nc")
dflow_map_3d = load_nc_info(datadir_3d, r"DCSM-FM_0_5nm_...._map.nc")

interp_2d = load_dflow_grid(dflow_map_2d, 50, true)
interp_3d = load_dflow_grid(dflow_map_3d, 50, true)

drifter = drifterdata(drifternumber)                                             # Retrieve data of the drifter
t0, starttime, endtime = driftertimes(drifter)

# Use the maximum end time which can be simulated
# Normally, the endtime should be the drifter endtime. However, drifter 2 endtime is later than the
# dflow map provides. If that is the case, the endtime of the dflow map should be used to avoid problems.
if water_3d
   endtime = min(endtime, get_times(dflow_map_3d,t0)[end])
else
   endtime = min(endtime, get_times(dflow_map_2d,t0)[end])
end


u_2d = initialize_interpolation(dflow_map_2d,interp_2d, "mesh2d_ucx",t0,0.0);
v_2d = initialize_interpolation(dflow_map_2d,interp_2d, "mesh2d_ucy",t0,0.0);
u_3d = initialize_interpolation(dflow_map_3d,interp_3d, "mesh2d_ucx",t0,0.0);
v_3d = initialize_interpolation(dflow_map_3d,interp_3d, "mesh2d_ucy",t0,0.0);
u_wind = initialize_interpolation(dflow_map_2d,interp_2d, "mesh2d_windx",t0,0.0);
v_wind = initialize_interpolation(dflow_map_2d,interp_2d, "mesh2d_windy",t0,0.0);
# In the 3D maps of March/April, wind data is also provided, thus then we load this one
# Else, we use the wind data from the 2D model (this should be equal)
if drifternumber <=4
   u_wind3 = initialize_interpolation(dflow_map_3d,interp_3d, "mesh2d_windx",t0,0.0)
   v_wind3 = initialize_interpolation(dflow_map_3d,interp_3d, "mesh2d_windy",t0,0.0)
end

wh=initialize_interpolation(matroos_data,"wave_height_hm0",t0,0.0)
wp=initialize_interpolation(matroos_data,"wave_period_tm10",t0,0.0)
wd=initialize_interpolation(matroos_data,"wave_dir_th0",t0,0.0)

variables=["lon","lat","age"]
d["variables"]=variables
m=length(variables)
p=zeros(m,n)
ds, s = track_of_drifter!(zeros(3), zeros(3), starttime, t0, 60, drifter)
p[1,:]=s[1]*ones(1,n)
p[2,:]=s[2]*ones(1,n)

d["particles"]=p #initial values
d["water_3d"] = water_3d

d["reftime"] = t0
d["dt"] = 300
d["tstart"] = starttime
d["tend"] = starttime+floor(Int64,(endtime-starttime)/d["dt"])*d["dt"]

###### Write to netcdf ######
d["write_maps_times"] = collect(starttime:3600:endtime)                          # Time at which data should be written to netcdf
d["write_maps"] = true
use2d3d = "3D"
if !water_3d
  use2d3d = "2D"
end
d["write_maps_filename"] = "netcdf_diffusie_drifter$(drifternumber)_$(use2d3d).nc"       # Save data in NetCDF file
d["write_maps_dir"]="case_drifters_german_bight/netcdf_output"

###### Plot maps ######
d["plot_maps_times"] = collect(starttime:14400:endtime)                          # Time at which plot should be made
d["plot_maps"] = true

d["plot_maps_folder"]="case_drifters_german_bight/images_drifter$(drifternumber)"

########################### Prepage background image ###########################
plot_maps_size = d["plot_maps_size"]
width,height = plot_maps_size
plot_bbox = d["bbox"]
wms_server = WmsServer("emodnet-bathymetry")
img = get_map(wms_server,plot_bbox,width,height)
d["background_image"] = img

function plot_background(d)
   plot_bbox = d["bbox"]
   img = d["background_image"]
   f = plot_image(img,plot_bbox)
   f = plot!(xaxis = ("Longitude \n ", (plot_bbox[1],plot_bbox[3]), font(30)), yaxis = (" \n Latitude", (plot_bbox[2],plot_bbox[4]), font(30)), legend = :topleft, legendfont = font(20), titlefontsize=25)
   f = plot!(convert(Array{Float64}, drifter[:,3]),convert(Array{Float64}, drifter[:,2]), linecolor = :black, legend=false)    # Plot track of drifter
   return(f)
end
d["plot_maps_background"] = plot_background

###### Velocity function for the particles ######
function f!(ds,s,t,i,d)
   x,y,age = s
   z=0.0
   R = 6371.0e3                                                                  # Mean radius of earth from wikipedia
   deg2rad = pi/180.0                                                            # Converts degrees to radians
   rad2deg = 180.0/pi                                                            # Converts radians to degrees
   up=0
   vp=0
   dt=d["dt"]
   uw=0
   vw=0
   ua=0
   va=0
   if water_3d
      uw = u_3d(x,y,z,t)
      vw = v_3d(x,y,z,t)
      if drifternumber <= 4
         ua = u_wind3(x,y,z,t)
         va = v_wind3(x,y,z,t)
      else
         ua = u_wind(x,y,z,t)
         va = v_wind(x,y,z,t)
      end
   else
      uw = u_2d(x,y,z,t)
      vw = v_2d(x,y,z,t)
      ua = u_wind(x,y,z,t)
      va = v_wind(x,y,z,t)
   end
   if i==1 #Use drifer data
      track_of_drifter!(ds,s,t,d["reftime"],d["dt"],drifter)
   elseif i<=301
      up = uw
      vp = vw
   elseif i<=601
      up = uw+0.016*ua*0.9
      vp = vw+0.016*va*0.9
   elseif i<=901
      usJ, vsJ = uv_sJ(wh(x,y,z,t),wp(x,y,z,t),wd(x,y,z,t))
      up = uw+usJ
      vp = vw+vsJ
   elseif i<=1201
      (up,vp) = water_stokes_wind(ua,va, uw,vw,0,0)
   elseif i<=1501
      (up,vp) = water_stokes_wind(ua,va,uw,vw,0.016*ua*0.5,0.016*va*0.5)
   else
      usJ, vsJ = uv_sJ(wh(x,y,z,t),wp(x,y,z,t),wd(x,y,z,t))
      (up,vp) = water_stokes_wind(ua,va, uw,vw,usJ,vsJ)
   end

   # Calculate and add turbulent diffusivity, using Pr=1
   if i>1
      # Estimate the Eddy viscosity and its derivates, using a Smagorinsky model
      if water_3d
         (K, Kdx, Kdy)=estimate_viscosity_smag(interp_3d, x,y,t,u_3d,v_3d)
      else
         (K, Kdx, Kdy)=estimate_viscosity_smag(interp_2d, x,y,t,u_2d,v_2d)
      end
      if !(uw==vw==ua==va==0.0)
         # https://doi.org/10.1016/j.ocemod.2017.11.008 eq. 27
         up += Kdy+randn()*sqrt(2*K*dt)/dt
         vp += Kdx+randn()*sqrt(2*K*dt)/dt
      end
   end

   # Convert the velocity in [m/s] to dlat and dlon for latitude and longitude
   ds[1] = dx   = rad2deg*up/(R*cos(deg2rad*s[2]))
   ds[2] = dy   = rad2deg*vp/R
   ds[3] = dage = 1.0
end
d["f"]=f!

run_simulation(d)
