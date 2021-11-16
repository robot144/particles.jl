#
# Run simulation for drifter model.
# Usage: julia --project case_drifter_challenge/drifters.jl /path/to/config.toml|csv"
#
# input should look something like this:
#id,x,y,nparticles,total_mass,radius,reftime,tend,dt,datapath
# 1,-69.93235,34.95782,10,5,0.02,2021-05-30T00:16:49,259200,0.0167,D:\particle_tracking_demo\FFT_Challenge\Model_validation\Data\CMEMS_data\hindcast_processed_withwind\

using Particles
using Plots

include(joinpath(@__DIR__, "drifterfunctions.jl"))

#
# read user input
#
usage = "Usage: julia --project=.. drifters.jl /path/to/config.toml|csv"
n = length(ARGS)
if n != 1
    throw(ArgumentError(usage))
end
config_path = only(ARGS)
if !isfile(config_path)
    throw(ArgumentError("File not found: $(config_path)\n" * usage))
else
    println("Reading config from file $(config_path).")
end
d = Particles.config(config_path)

#
# check and iniitialize
#
if !isa(d["id"], Vector)
    d["id"] = [d["id"]]
    d["x"] = [d["x"]]
    d["y"] = [d["y"]]
end
d["nsources"] = length(d["id"])
if length(d["x"]) != d["nsources"]
    error("Length of x not equal to $(d["nsources"])")
end
if length(d["y"]) != d["nsources"]
    error("Length of y not equal to $(d["nsources"])")
end

if !haskey(d, "coordinates")
    d["coordinates"] = "spherical"
end

d["time_direction"] = :forwards # :forwards or :backwards

#
# flow data from delft3d-fm
# optionally also winds through this route
#
#dflow_map = load_nc_info(d["datapath"], r"DCSM-FM_05nm_...._map.nc")
#const interp = load_dflow_grid(dflow_map, 50, true)

#u = initialize_interpolation(dflow_map, interp, "mesh2d_ucx", d["reftime"], 0.0, d["time_direction"]);
#v = initialize_interpolation(dflow_map, interp, "mesh2d_ucy", d["reftime"], 0.0, d["time_direction"]);
#u_wind = initialize_interpolation(dflow_map, interp, "mesh2d_windx", d["reftime"], 0.0, d["time_direction"])
#v_wind = initialize_interpolation(dflow_map, interp, "mesh2d_windy", d["reftime"], 0.0, d["time_direction"])

# flow data from cmems
cmems_u = CmemsData(d["datapath"], "u0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc")
cmems_v = CmemsData(d["datapath"], "v0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc")
t0 = d["reftime"]
u = initialize_interpolation(cmems_u, "uo", t0, NaN)  # water velocity x-dir
v = initialize_interpolation(cmems_v, "vo", t0, NaN)  # water velocity y-dir

# wind data from gfs
gfs_u = GFSData(d["datapath"], "gfs_winds.nc")
gfs_v = GFSData(d["datapath"], "gfs_winds.nc")
t0 = d["reftime"]
u_wind = initialize_interpolation(gfs_u, "10u", t0, NaN)  # wind velocity x-dir
v_wind = initialize_interpolation(gfs_v, "10v", t0, NaN)  # wind velocity y-dir



variables = ["lon", "lat", "age"]
d["variables"] = variables
m = length(variables)
n = d["nparticles"] * d["nsources"]
p = zeros(m, n)
# ds, s = track_of_drifter!(zeros(3), zeros(3), starttime, t0, 60, drifter)
#p[1,:] .= rand(-1.:0.00001:1., d["nparticles"]) .* d["radius"] .+ d["x"]
#p[2,:] .= rand(-1.:0.00001:1., d["nparticles"]) .* d["radius"] .+ d["y"]
iindex = 1
for isource = 1:d["nsources"]
    xsource = d["x"][isource]
    ysource = d["y"][isource]
    dx = d["radius"]
    dy = dx
    for ipart = 1:d["nparticles"]
        global iindex
        p[1, iindex] = xsource + dx * randn()
        p[2, iindex] = ysource + dy * randn()
        iindex += 1
    end
end
d["ids"] = repeat(d["id"], inner = d["nparticles"])  # repeated ids for each particle
d["nparticles"] = n #NOTE overwrite from input

d["particles"] = p # initial values

@info p

###### Write to netcdf ######
haskey(d, "write_maps") || (d["write_maps"] = true)
haskey(d, "write_maps_filename") || (d["write_maps_filename"] = "drifters.nc")    # Save data in NetCDF file
haskey(d, "write_maps_dir") || (d["write_maps_dir"] = "netcdf_output")
d["write_maps_times"] = collect(d["tstart"]:3600:d["tend"])                          # Time at which data should be written to netcdf

###### Plot maps ######
haskey(d, "plot_maps") || (d["plot_maps"] = true)
haskey(d, "plot_maps_folder") || (d["plot_maps_folder"] = "images_drifter")
d["plot_maps_times"] = collect(d["tstart"]:(24*3600):d["tend"])         # Time at which plot should be made
# spatial domain for plots
if !haskey(d, "bbox")
    dx = maximum(d["x"]) - minimum(d["x"])
    dy = maximum(d["y"]) - minimum(d["y"])
    d["buffer"] = d["radius"] * 100
    d["bbox"] = [minimum(d["x"]) - d["buffer"], minimum(d["y"]) - d["buffer"], maximum(d["x"]) + d["buffer"], maximum(d["y"]) + d["buffer"]]                                                 # Where we expect particles
    d["plot_maps_size"] = (1500, 1500)  # base this on bbox
end


@info "Using the following  " d

########################### Prepare background image ###########################
plot_maps_size = d["plot_maps_size"]
width, height = plot_maps_size
plot_bbox = d["bbox"]
if haskey(d, "plot_background_source")
    wms_server = WmsServer(d["plot_background_source"])
else
    #wms_server = WmsServer("emodnet-bathymetry")
    wms_server = WmsServer("gebco")
end
img = get_map(wms_server, plot_bbox, width, height)
d["background_image"] = img

function plot_background(d)
    plot_bbox = d["bbox"]
    img = d["background_image"]
    f = plot_image(img, plot_bbox)
    f = plot!(xaxis = ("Longitude \n ", (plot_bbox[1], plot_bbox[3]), font(30)), yaxis = (" \n Latitude", (plot_bbox[2], plot_bbox[4]), font(30)), legend = :topleft, legendfont = font(20), titlefontsize = 25)
    return (f)
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
    # up = uw
    # vp = vw
    # 2: Flow plus a factor times wind
    up = uw + ua * d["leeway_coeff"]
    vp = vw + va * d["leeway_coeff"]
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
    #(K, Kdx, Kdy) = estimate_viscosity_smag(interp, x, y, t, u, v)
    # WORKAROUND
    K = d["K"]
    Kdx = d["Kdx"]
    Kdy = d["Kdy"]
    # TODO: fix estimate_viscosity_smag for regular grids
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
