#
# Run simulation for drifter model.
# Usage: julia --project case_drifter_challenge/drifters.jl /path/to/config.toml|csv"
#

cd("../Particles.jl")
import Pkg
using Pkg

Pkg.activate(".")
Pkg.instantiate()

using Particles
using Plots

#include(joinpath(@__DIR__, "drifterfunctions.jl"))

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
# check and initialize
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

# simulation timing
d["time_direction"] = :forwards # :forwards or :backwards
if !haskey(d, "start")
    d["tstart"] = 0.0
end
if !haskey(d, "tend")
    error("Final time of simulation tend is missing in config.")
end
if !haskey(d, "dt")
    error("Time-step in seconds dt is missing in config.")
end
if !haskey(d, "write_maps_interval")
    d["write_maps_interval"] = 1800.0 # seconds
end
if !haskey(d, "plot_maps_interval")
    d["plot_maps_interval"] = 7200.0 # seconds
end

#
# flow data from delft3d-fm
# optionally also winds through this route
#

function zero_fun(x, y, z, t) #zero everywhere
    return 0.0
end


# check keywords for flowdata
if !haskey(d, "current_dir")
    d["current_dir"] = d["datapath"]
end
current_dir = d["current_dir"]
if haskey(d, "current_filename")
    d["current_x_filename"] = d["current_filename"]
    d["current_y_filename"] = d["current_filename"]
end
if !haskey(d, "current_x_filename")
    error("Missing key: current_x_filename")
end
if !haskey(d, "current_y_filename")
    error("Missing key: current_y_filename")
end
if !haskey(d, "current_filetype")
    error("Missing key: current_filetype")
end
if !haskey(d,"current_x_var")
    d["current_x_var"] = "longitude"
end
if !haskey(d,"current_y_var")
    d["current_y_var"] = "latitude"
end
if lowercase(d["current_filetype"]) == "cmems"
    if !haskey(d,"current_ucx_var")
        d["current_ucx_var"] = "uo"
    end
    if !haskey(d,"current_ucy_var")
        d["current_ucy_var"] = "vo"
    end
elseif lowercase(d["current_filetype"]) == "delft3d-fm"
    if !haskey(d,"current_ucx_var")
        d["current_ucx_var"] = "mesh2d_ucx"
    end
    if !haskey(d,"current_ucy_var")
        d["current_ucy_var"] = "mesh2d_ucy"
    end
end

# create u and v functions for flowdata
if lowercase(d["current_filetype"]) == "cmems"
    cmems_u = CmemsData(current_dir, d["current_x_filename"]; lon = d["current_x_var"], lat = d["current_y_var"])
    cmems_v = CmemsData(current_dir, d["current_y_filename"]; lon = d["current_x_var"], lat = d["current_y_var"])
    t0 = d["reftime"]
    u = initialize_interpolation(cmems_u, d["current_ucx_var"], t0, 0.0)  # water velocity x-dir
    v = initialize_interpolation(cmems_v, d["current_ucy_var"], t0, 0.0)  # water velocity y-dir
elseif lowercase(d["current_filetype"]) == "delft3d-fm"
    dflow_map = load_nc_info(current_dir, d["current_filename"])
    const interp = load_dflow_grid(dflow_map, 50, true)
    u = initialize_interpolation(dflow_map, interp, d["current_ucx_var"], d["reftime"], 0.0, d["time_direction"]);
    v = initialize_interpolation(dflow_map, interp, d["current_ucy_var"], d["reftime"], 0.0, d["time_direction"]);
elseif lowercase(d["current_filetype"]) == "zero"
    u = zero_fun
    v = zero_fun
else
    error("Invalid current_filtype: $(d["current_filetype"])")
end

# create u and v functions for SECOND flowdata
if haskey(d, "current2_filetype") 
    if lowercase(d["current2_filetype"]) == "delft3d-fm"
        dflow_map = load_nc_info(d["current2_dir"], d["current2_filename"])
        const interp = load_dflow_grid(dflow_map, 50, true)
        u2 = initialize_interpolation(dflow_map, interp, d["current2_ucx_var"], d["reftime"], 0.0, d["time_direction"]);
        v2 = initialize_interpolation(dflow_map, interp, d["current2_ucy_var"], d["reftime"], 0.0, d["time_direction"]);
        println("A second flow-field will be used from: $(d["current2_filename"])")
    elseif lowercase(d["current2_filetype"]) == "zero"
        u2 = zero_fun
        v2 = zero_fun
    else
        error("Invalid current2_filetype (only 'delft3d-fm' is supported): $(d["current2_filetype"])")
    end
else
    u2 = zero_fun
    v2 = zero_fun
end

# check input for winds
if !haskey(d, "wind_dir")
    d["wind_dir"] = d["datapath"]
end
if haskey(d, "wind_filename")
    d["wind_x_filename"] = d["wind_filename"]
    d["wind_y_filename"] = d["wind_filename"]
end
if !haskey(d, "wind_x_filename")
    error("Missing key: wind_x_filename")
end
if !haskey(d, "wind_y_filename")
    error("Missing key: wind_y_filename")
end
if !haskey(d, "wind_filetype")
    error("Missing key: wind_filetype")
end
if !haskey(d,"wind_x_var")
   d["wind_x_var"] = "x"
end
if !haskey(d,"wind_y_var")
   d["wind_y_var"] = "y"
end

# create u_wind and v_wind
if lowercase(d["wind_filetype"]) == "gfs"
    # wind data from gfs
    gfs_u = GFSData(d["wind_dir"], d["wind_x_filename"]; lon = d["wind_x_var"], lat = d["wind_y_var"])
    gfs_v = GFSData(d["wind_dir"], d["wind_y_filename"]; lon = d["wind_x_var"], lat = d["wind_y_var"])
    t0 = d["reftime"]
    if !haskey(d,"wind_x_wrap") 
        if all(d["x"] .>= -180.0 * d["x"] .<= 180.0) && all(d["y"] .>= -90.0 * d["y"] .<= 90.0) && all(gfs_u.grid.xnodes .>= 0.0 * gfs_u.grid.xnodes .<= 360.0)
            @warn "The particles seem to be in spherical coordinates (lon: [-180 180]), but the provided GFS data is in lon: [0 360]. wind_x_wrap is set to 'true'"
            d["wind_x_wrap"] = true  # shift longitudes from [0 360] to [-180 180]
        else
            d["wind_x_wrap"] = false
        end
    end
    u_wind = initialize_interpolation(gfs_u, "10u", t0, NaN, wrap = d["wind_x_wrap"])  # wind velocity x-dir
    v_wind = initialize_interpolation(gfs_v, "10v", t0, NaN, wrap = d["wind_x_wrap"])  # wind velocity y-dir
elseif lowercase(d["wind_filetype"]) == "delft3d-fm"
    u_wind = initialize_interpolation(dflow_map, interp, "mesh2d_windx", d["reftime"], 0.0, d["time_direction"])
    v_wind = initialize_interpolation(dflow_map, interp, "mesh2d_windy", d["reftime"], 0.0, d["time_direction"])
    error("TODO: make this work.")
elseif lowercase(d["wind_filetype"]) == "zero"
    @warn "Running without wind"
    u_wind = zero_fun
    v_wind = zero_fun
else
    error("Invalid wind_filtype: $(d["current_filetype"])")
end


# initialize
variables = ["lon", "lat", "age"]
d["variables"] = variables
m = length(variables)
n = d["npartpersource"] * d["nsources"]
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
    for ipart = 1:d["npartpersource"]
        global iindex
        p[1, iindex] = xsource + dx * randn()
        p[2, iindex] = ysource + dy * randn()
        iindex += 1
    end
end
d["ids"] = repeat(d["id"], inner = d["npartpersource"])  # repeated ids for each particle
d["nparticles"] = n #NOTE total nr of particles (=nsources*npartpersource)

d["particles"] = p # initial values

@info p

###### Write to netcdf ######
haskey(d, "write_maps") || (d["write_maps"] = true)
if d["write_maps"]
    haskey(d, "write_maps_filename") || (d["write_maps_filename"] = "drifters.nc")    # Save data in NetCDF file
    haskey(d, "write_maps_dir") || (d["write_maps_dir"] = "netcdf_output")
    d["write_maps_times"] = collect(d["tstart"]:d["write_maps_interval"]:d["tend"])                          # Time at which data should be written to netcdf
end

###### Plot maps ######
haskey(d, "plot_maps") || (d["plot_maps"] = true)
if d["plot_maps"]
    haskey(d, "plot_maps_folder") || (d["plot_maps_folder"] = "images_drifter")
    d["plot_maps_times"] = collect(d["tstart"]:(d["plot_maps_interval"]):d["tend"])         # Time at which plot should be made
    # spatial domain for plots
    if !haskey(d, "bbox")
        dx = maximum(d["x"]) - minimum(d["x"])
        dy = maximum(d["y"]) - minimum(d["y"])
        d["buffer"] = d["radius"] * 100
        d["bbox"] = [minimum(d["x"]) - d["buffer"], minimum(d["y"]) - d["buffer"], maximum(d["x"]) + d["buffer"], maximum(d["y"]) + d["buffer"]]                                                 # Where we expect particles
        d["plot_maps_size"] = (1500, 1500)  # base this on bbox
    end
end

@info "Using the following  " d

########################### Prepare background image ###########################
if d["plot_maps"]
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
end

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
    uw += u2(x, y, z, t)
    vw += v2(x, y, z, t) 
    if uw != 0
        ua = u_wind(x, y, z, t)
        va = v_wind(x, y, z, t)
        up = uw + ua * d["leeway_coeff"]
        vp = vw + va * d["leeway_coeff"]
    end

    # Various models:
    # 0: Use drifer data
    # track_of_drifter!(ds,s,t,d["reftime"],d["dt"],drifter)
    # 1: Only flow velocities
    # up = uw
    # vp = vw
    # 2: Flow plus a factor times wind
    # up = uw + ua * d["leeway_coeff"]
    # vp = vw + va * d["leeway_coeff"]
    # 3: Add stokes drift to flow velocity
    # usJ, vsJ = uv_sJ(wh(x,y,z,t),wp(x,y,z,t),wd(x,y,z,t))
    # up = uw+usJ
    # vp = vw+vsJ
    # 4: Combine flow, stokes drift and wind in an equilibrium for the particle velocity
    # usJ, vsJ = uv_sJ(wh(x,y,z,t),wp(x,y,z,t),wd(x,y,z,t))
    # (up,vp) = water_stokes_wind(ua,va, uw,vw,usJ,vsJ)
    # 5: Flow plus flow plus a factor times wind (a second flow field to allow for e.g. CMEMS + GTSM flow)

    # Calculate and add turbulent diffusivity, using Pr=1
    # Estimate the Eddy viscosity and its derivates, using a Smagorinsky model
    # This is only because the horizontal diffusion is not in the flow output files
    #(K, Kdx, Kdy) = estimate_viscosity_smag(interp, x, y, t, u, v)
    # WORKAROUND
    K = d["K"]
    Kdx = 0.0
    Kdy = 0.0
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

if d["write_maps"] && haskey(d, "npartpersource") && d["npartpersource"] > 1
    import NetCDF
    using NetCDF
    import Statistics
    using Statistics
    nsources = d["nsources"]
    npartpersource = d["npartpersource"]

    fullfile = joinpath(d["write_maps_dir"], d["write_maps_filename"])
    file = NetCDF.open(fullfile)
    gatts = file.gatts
    time = ncread(fullfile, "time")
    ntimes = length(time)
    time_atts = file["time"].atts
    finalize(file)
    
    fullfile_mean = joinpath(d["write_maps_dir"], "source-averaged_" * d["write_maps_filename"])
    if isfile(fullfile_mean)
        println("Source-averaged output file exists. Removing file $(fullfile_mean)")
        rm(fullfile_mean)
    end
    nc = NetCDF.create(fullfile_mean, gatts = gatts, mode = NC_NETCDF4)

    for varname in keys(file.vars)
        dimnames = [file[varname].dim[i].name for i = 1:file[varname].ndim]
        if "time" in dimnames && "particles" in dimnames
            data = ncread(fullfile, varname)
            nccreate(fullfile_mean, varname, "time", time, "sources", collect(1:1:nsources))
            data_mean = zeros(ntimes, nsources)
            for srci = 1:nsources
                ind1 = (srci - 1) * npartpersource + 1
                ind2 = srci * npartpersource
                data_mean[:,srci] = mean(data[:,ind1:ind2], dims = 2) #todo: take nanmean or skipmissing to avoid mean([1 2 NaN/missing 5]) to become NaN
            end
            ncwrite(data_mean, fullfile_mean, varname)
            ncputatt(fullfile_mean, varname, file[varname].atts)
        end
    end
    ncputatt(fullfile_mean, "time", time_atts)
    finalize(fullfile_mean)
end
