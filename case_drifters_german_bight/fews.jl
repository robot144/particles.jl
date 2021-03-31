using Particles
using Plots

include(joinpath(@__DIR__, "drifterfunctions.jl"))

usage = "Usage: julia --project case_drifters_german_bight/fews.jl /path/to/config.toml|csv"
n = length(ARGS)
if n != 1
    throw(ArgumentError(usage))
end
config_path = only(ARGS)
if !isfile(config_path)
    throw(ArgumentError("File not found: $(config_path)\n" * usage))
end
d = Particles.config(config_path)

d["coordinates"] = "spherical"
d["buffer"] = d["radius"] * 100
d["bbox"] = [d["x"] - d["buffer"],d["y"] - d["buffer"],d["x"] + d["buffer"],d["y"] + d["buffer"]]                                                 # Where we expect particles
d["plot_maps_size"] = (1500, 1500)  # base this on bbox
d["time_direction"] = :forwards # :forwards or :backwards

dflow_map = load_nc_info(d["datapath"], r"DCSM-FM_05nm_...._map.nc")
const interp = load_dflow_grid(dflow_map, 50, true)

u = initialize_interpolation(dflow_map, interp, "mesh2d_ucx", d["reftime"], 0.0, d["time_direction"]);
v = initialize_interpolation(dflow_map, interp, "mesh2d_ucy", d["reftime"], 0.0, d["time_direction"]);
u_wind = initialize_interpolation(dflow_map, interp, "mesh2d_windx", d["reftime"], 0.0, d["time_direction"])
v_wind = initialize_interpolation(dflow_map, interp, "mesh2d_windy", d["reftime"], 0.0, d["time_direction"])

variables = ["lon","lat","age"]
d["variables"] = variables
m = length(variables)
p = zeros(m, d["nparticles"])
# ds, s = track_of_drifter!(zeros(3), zeros(3), starttime, t0, 60, drifter)
p[1,:] .= rand(-1.:0.00001:1., d["nparticles"]) .* d["radius"] .+ d["x"]
p[2,:] .= rand(-1.:0.00001:1., d["nparticles"]) .* d["radius"] .+ d["y"]
d["particles"] = p # initial values
@info p

###### Write to netcdf ######
d["write_maps_times"] = collect(d["tstart"]:3600:d["tend"])                          # Time at which data should be written to netcdf
d["write_maps"] = true
d["write_maps_filename"] = "netcdf_diffusie_drifter$(d["id"]).nc"       # Save data in NetCDF file
d["write_maps_dir"] = "netcdf_output"

###### Plot maps ######
d["plot_maps_times"] = collect(d["tstart"]:(24 * 3600):d["tend"])                          # Time at which plot should be made
d["plot_maps"] = true

d["plot_maps_folder"] = "images_drifter"
@info "Using the following config: " d

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
   ds.x = rad2deg * up / (R * cos(deg2rad * s[2]))
   ds.y = rad2deg * vp / R
   ds.z = 1.0
end
d["f"] = f!

run_simulation(d)
