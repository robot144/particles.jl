# Interact with dflow netcdf output map files
#
using Dates
using Zarr

#
# test
#

function test_estuary() #2D maps of a simple estuary. This model has only two domains
   #initialize
   if isfile(joinpath(pwd(),"config_maps_interp.toml"))
      rm(joinpath(pwd(),"config_maps_interp.toml"),recursive=false)
   end
   if isdir(joinpath(pwd(),"estuary_map.zarr"))
      rm(joinpath(pwd(),"estuary_map.zarr"),recursive=true)
   end

   # generate config file
   netcdf_filenames = ["../test_data/estuary_0000_map.nc", "../test_data/estuary_0001_map.nc"]
   include("../src/dflow_map_interp_to_zarr.jl")
   #main(netcdf_filenames)
   result=Base.invokelatest(main,netcdf_filenames) #invokelatest is needed for include in function

   configfile=joinpath(pwd(),"config_maps_interp.toml")
   @test isfile(configfile)
   # load configfile
   config=TOML.parsefile(configfile)
   @test haskey(config,"global")
   @test haskey(config,"waterlevel")
   @test haskey(config,"x_velocity")

   # run interpolation with default config file
   result=Base.invokelatest(main,[configfile]) #invokelatest is needed for include in function

   # check results
   output_folder=joinpath(pwd(),"estuary_map.zarr")
   @test isdir(output_folder)
   z = zopen(output_folder)

   # check variables
   @test haskey(z,"x_center")
   xs = z["x_center"][:]
   @test length(xs) == 200
   @test haskey(z,"y_center")
   ys = z["y_center"][:]
   @test length(ys) == 1
   @test haskey(z,"x_velocity")
   u = z["x_velocity"][:,:,:]
   @test size(u) == (200,1,25)
   @test haskey(z,"y_velocity")
   v = z["y_velocity"][:,:,:]
   @test size(v) == (200,1,25)
   @test haskey(z,"waterlevel")
   h = z["waterlevel"][:,:,:]
   @test size(h) == (200,1,25)
   @test haskey(z,"time")
   t = z["time"]
   @test length(t[:]) == 25
   @test t.attrs["units"] == "seconds since 1991-01-01 00:00:00"
   @test t.attrs["add_offset"] ≈ 0.0
   @test t.attrs["scale_factor"] ≈ 1.0
   @test t.attrs["standard_name"] == "time" 


   #finalize
   # if isfile(joinpath(pwd(),"config_maps_interp.toml"))
   #    rm(joinpath(pwd(),"config_maps_interp.toml"),recursive=false)
   # end
   # if isdir(joinpath(pwd(),"estuary_map.zarr"))
   #    rm(joinpath(pwd(),"estuary_map.zarr"),recursive=true)
   # end
end

function test_f34() #2D maps of a tidal inlet. This model has only one domain, but has some real geometry and dry cells
   #initialize
   if isfile(joinpath(pwd(),"config_maps_interp.toml"))
      rm(joinpath(pwd(),"config_maps_interp.toml"),recursive=false)
   end
   if isdir(joinpath(pwd(),"f34_map.zarr"))
      rm(joinpath(pwd(),"f34_map.zarr"),recursive=true)
   end

   # generate config file
   netcdf_filenames = ["../test_data/f34_map.nc"]
   include("../src/dflow_map_interp_to_zarr.jl")
   #main(netcdf_filenames)
   result=Base.invokelatest(main,netcdf_filenames) #invokelatest is needed for include in function

   configfile=joinpath(pwd(),"config_maps_interp.toml")
   @test isfile(configfile)
   # load configfile
   config=TOML.parsefile(configfile)
   @test haskey(config,"global")
   @test haskey(config,"waterlevel")
   @test haskey(config,"x_velocity")

   # run interpolation with default config file
   result=Base.invokelatest(main,[configfile]) #invokelatest is needed for include in function

   # check results
   output_folder=joinpath(pwd(),"f34_map.zarr")
   @test isdir(output_folder)
   z = zopen(output_folder)

   # check variables
   @test haskey(z,"x_center")
   xs = z["x_center"][:]
   @test length(xs) == 14
   @test haskey(z,"y_center")
   ys = z["y_center"][:]
   @test length(ys) == 11
   @test haskey(z,"x_velocity")
   u = z["x_velocity"][:,:,:]
   @test size(u) == (14,11,151)
   @test haskey(z,"y_velocity")
   v = z["y_velocity"][:,:,:]
   @test size(v) == (14,11,151)
   @test haskey(z,"waterlevel")
   h = z["waterlevel"][:,:,:]
   @test size(h) == (14,11,151)
   h12=h[:,:,12] #check dry cells
   @test h12[1,1] == 9999
   println("h=$(h12)")
   @test haskey(z,"time")
   t = z["time"]
   @test length(t[:]) == 151
   @test t.attrs["units"] == "seconds since 1990-08-05 00:00:00 +00:00"
   @test t.attrs["add_offset"] ≈ 0.0
   @test t.attrs["scale_factor"] ≈ 1.0
   @test t.attrs["standard_name"] == "time" 


   #finalize
   # if isfile(joinpath(pwd(),"config_maps_interp.toml"))
   #    rm(joinpath(pwd(),"config_maps_interp.toml"),recursive=false)
   # end
   # if isdir(joinpath(pwd(),"estuary_map.zarr"))
   #    rm(joinpath(pwd(),"estuary_map.zarr"),recursive=true)
   # end
end

function test_lock_exchange() 
   #3D maps of a lock exchange with z-layers. This model has only 1 domain
   # The netcdf file contains fullgrid arrays for the layers
   #initialize
   if isfile(joinpath(pwd(),"config_maps_interp.toml"))
      rm(joinpath(pwd(),"config_maps_interp.toml"),recursive=false)
   end
   if isdir(joinpath(pwd(),"locxxz_fullgrid_map.zarr"))
      rm(joinpath(pwd(),"locxxz_fullgrid_map.zarr"),recursive=true)
   end

   # generate config file
   netcdf_filenames = ["../test_data/locxxz_fullgrid_map.nc"]
   include("../src/dflow_map_interp_to_zarr.jl")
   #main(netcdf_filenames)
   result=Base.invokelatest(main,netcdf_filenames) #invokelatest is needed for include in function

   configfile=joinpath(pwd(),"config_maps_interp.toml")
   @test isfile(configfile)
   # load configfile
   config=TOML.parsefile(configfile)
   @test haskey(config,"global")
   @test haskey(config,"waterlevel")
   @test haskey(config,"x_velocity")

   @test haskey(config,"z_velocity")

   # run interpolation with default config file
   result=Base.invokelatest(main,[configfile]) #invokelatest is needed for include in function

   # check results
   output_folder=joinpath(pwd(),"locxxz_fullgrid_map.zarr")
   @test isdir(output_folder)
   z = zopen(output_folder)

   # check variables
   @test haskey(z,"x_center")
   xs = z["x_center"][:]
   @test length(xs) == 344
   @test haskey(z,"y_center")
   ys = z["y_center"][:]
   @test length(ys) == 1
   @test haskey(z,"x_velocity")
   u = z["x_velocity"][:,:,:,:] # x y z t
   @test size(u) == (344,1,40,16)
   @test haskey(z,"y_velocity")
   v = z["y_velocity"][:,:,:,:]
   @test size(v) == (344,1,40,16)
   @test haskey(z,"z_velocity")
   w = z["z_velocity"][:,:,:,:]
   @test size(w) == (344,1,41,16) # one more interface than layers

   @test haskey(z,"eddy_visc_z")
   nu = z["eddy_visc_z"][:,:,:,:]
   @test size(nu) == (344,1,41,16) # one more interface than layers  

   @test haskey(z,"waterlevel")
   h = z["waterlevel"][:,:,:]
   @test size(h) == (344,1,16)
   @test haskey(z,"time")
   t = z["time"]
   @test length(t[:]) == 16
   @test t.attrs["units"] == "seconds since 2001-01-01 00:00:00 +00:00"
   @test t.attrs["add_offset"] ≈ 0.0
   @test t.attrs["scale_factor"] ≈ 1.0
   @test t.attrs["standard_name"] == "time" 
   # "z_center_3d" and  "z_iface_3d" 
   @test haskey(z,"z_center_3d")
   zc = z["z_center_3d"][:,:,:,:] # x y z t
   @test size(zc) == (344,1,40,16)
   @test haskey(z,"z_iface_3d")
   zw = z["z_iface_3d"][:,:,:,:] # x y z t
   @test size(zw) == (344,1,41,16)


   #finalize
   # if isfile(joinpath(pwd(),"config_maps_interp.toml"))
   #    rm(joinpath(pwd(),"config_maps_interp.toml"),recursive=false)
   # end
   # if isdir(joinpath(pwd(),"estuary_map.zarr"))
   #    rm(joinpath(pwd(),"estuary_map.zarr"),recursive=true)
   # end
end


test_estuary()
test_f34()
test_lock_exchange()

