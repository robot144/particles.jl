using Particles
using Plots
include("drifterfunctions.jl")

drifternumber = 4;

d=default_userdata()
n=10 #number of particles

d["nparticles"]=n
d["coordinates"]="spherical"
d["bbox"] = [6.0,53.5,9.1,55.5]                                                 # Where we expect particles
d["plot_maps_size"] = (1500, 1000)
d["time_direction"] = :forwards # :forwards or :backwards


# Use different data directories, depending on the drifternumber and water_3d
# Drifter 1-4 are in the period of march/april of 2017, thus dflow and matroos data of this period should be used
# Drifter 5-7 are in the period of october of 2017. The data for this run is currently not available here.
if drifternumber >= 5
   datadir = "data/dcsm-fm_201710"
else
   datadir = "data/dcsm-fm_201703"
end

dflow_map = load_nc_info(datadir, r"DCSM-FM_0_5nm_...._map.nc")
interp = load_dflow_grid(dflow_map, 50, true)

drifter = drifterdata(drifternumber)                                             # Retrieve data of the drifter
t0, starttime, endtime = driftertimes(drifter)

# Use the maximum end time which can be simulated
# Normally, the endtime should be the drifter endtime. However, drifter 2 endtime is later than the
# dflow map provides. If that is the case, the endtime of the dflow map should be used to avoid problems.
endtime = min(endtime, get_times(dflow_map, t0)[end])

u = initialize_interpolation(dflow_map, interp, "mesh2d_ucx", t0, 0.0, d["time_direction"]);
v = initialize_interpolation(dflow_map, interp, "mesh2d_ucy", t0, 0.0, d["time_direction"]);

# u_wind = initialize_interpolation(dflow_map,interp_2d, "mesh2d_windx",t0,0.0,d["time_direction"]);
# v_wind = initialize_interpolation(dflow_map,interp_2d, "mesh2d_windy",t0,0.0,d["time_direction"]);
# In the 3D maps of March/April, wind data is also provided, thus then we load this one
# Else, we use the wind data from the 2D model (this should be equal)
if drifternumber <= 4
   u_wind = initialize_interpolation(dflow_map, interp, "mesh2d_windx", t0, 0.0, d["time_direction"])
   v_wind = initialize_interpolation(dflow_map, interp, "mesh2d_windy", t0, 0.0, d["time_direction"])
end

# wave data
# wh=initialize_interpolation(matroos_data,"wave_height_hm0",t0,0.0,d["time_direction"])
# wp=initialize_interpolation(matroos_data,"wave_period_tm10",t0,0.0,d["time_direction"])
# wd=initialize_interpolation(matroos_data,"wave_dir_th0",t0,0.0,d["time_direction"])

variables = ["lon","lat","age"]
d["variables"] = variables
m = length(variables)
p = zeros(m, n)
ds, s = track_of_drifter!(zeros(3), zeros(3), starttime, t0, 60, drifter)
p[1,:] = s[1] * ones(1, n)
p[2,:] = s[2] * ones(1, n)

d["particles"] = p # initial values
d["reftime"] = t0
d["dt"] = 300
d["tstart"] = starttime
d["tend"] = endtime

###### Write to netcdf ######
d["write_maps_times"] = collect(starttime:3600:endtime)                          # Time at which data should be written to netcdf
d["write_maps"] = true
use2d3d = "3D"
d["write_maps_filename"] = "netcdf_diffusie_drifter$(drifternumber)_$(use2d3d).nc"       # Save data in NetCDF file
d["write_maps_dir"] = "netcdf_output"

###### Plot maps ######
d["plot_maps_times"] = collect(starttime:(24 * 3600):endtime)                          # Time at which plot should be made
d["plot_maps"] = true

d["plot_maps_folder"] = "images_drifter$(drifternumber)"

########################### Prepare background image ###########################
plot_maps_size = d["plot_maps_size"]
width, height = plot_maps_size
plot_bbox = d["bbox"]
wms_server = WmsServer("emodnet-bathymetry")
img = get_map(wms_server, plot_bbox, width, height)
d["background_image"] = img

function plot_background(d)
   plot_bbox = d["bbox"]
   img = d["background_image"]
   f = plot_image(img, plot_bbox)
   f = plot!(xaxis=("Longitude \n ", (plot_bbox[1], plot_bbox[3]), font(30)), yaxis=(" \n Latitude", (plot_bbox[2], plot_bbox[4]), font(30)), legend=:topleft, legendfont=font(20), titlefontsize=25)
   f = plot!(convert(Array{Float64}, drifter[:,3]), convert(Array{Float64}, drifter[:,2]), linecolor=:black, legend=false)    # Plot track of drifter
   return(f)
end
d["plot_maps_background"] = plot_background

###### Velocity function for the particles ######
function f!(ds, s, t, i, d)
   x, y, age = s
   z = 0.0
   R = 6371.0e3                                                                  # Mean radius of earth from wikipedia
   deg2rad = pi / 180.0                                                            # Converts degrees to radians
   rad2deg = 180.0 / pi                                                            # Converts radians to degrees
   up = 0
   vp = 0
   dt = d["dt"]
   uw = 0
   vw = 0
   ua = 0
   va = 0
   uw = u(x, y, z, t)
   vw = v(x, y, z, t)
   ua = u_wind(x, y, z, t)
   va = v_wind(x, y, z, t)

   # Various models:
   # 0: Use drifer data
   # track_of_drifter!(ds,s,t,d["reftime"],d["dt"],drifter)
   # 1: Only flow velocities
   up = uw
   vp = vw
   # 2: Flow plus a factor times wind
   # up = uw+0.016*ua*0.9
   # vp = vw+0.016*va*0.9
   # 3: Add stokes drift to flow velocity
   # usJ, vsJ = uv_sJ(wh(x,y,z,t),wp(x,y,z,t),wd(x,y,z,t))
   # up = uw+usJ
   # vp = vw+vsJ
   # 4: Combine flow, stokes drift and wind in an equilibrium for the particle velocity
   # usJ, vsJ = uv_sJ(wh(x,y,z,t),wp(x,y,z,t),wd(x,y,z,t))
   # (up,vp) = water_stokes_wind(ua,va, uw,vw,usJ,vsJ)

   # Calculate and add turbulent diffusivity, using Pr=1
   # Estimate the Eddy viscosity and its derivates, using a Smagorinsky model
   # This is only because the horizontal diffusion is not in the flow output files
   (K, Kdx, Kdy) = estimate_viscosity_smag(interp, x, y, t, u, v)
   if !(uw == vw == ua == va == 0.0)
      # https://doi.org/10.1016/j.ocemod.2017.11.008 eq. 27
      up += Kdy + randn() * sqrt(2 * K * dt) / dt
      vp += Kdx + randn() * sqrt(2 * K * dt) / dt
   end
   if d["time_direction"] == :backwards
      up *= -1
      vp *= -1
   end

   # Convert the velocity in [m/s] to dlat and dlon for latitude and longitude
   ds[1] = dx   = rad2deg*up/(R*cos(deg2rad*s[2]))
   ds[2] = dy   = rad2deg*vp/R
   ds[3] = dage = 1.0
end
d["f"]=f!

run_simulation(d)
