using Test, Particles, Printf, NetCDF

#clear cache
if isdir(joinpath(pwd(),".cache"))
   rm(joinpath(pwd(),".cache"),recursive=true)
end

@testset "Particle tests" begin
    
   @testset "1d interpolation" begin
      include("test_1d_interpolation.jl")
   end

   @testset "Grid indexing" begin
      include("test_unstructured_grid.jl")
   end

   @testset "Cartesian grids " begin
      include("test_cartesian_grid.jl")
   end
      
   @testset "Background images over WMS" begin
      include("test_wms_client.jl")
   end

   @testset "Read and interpolate map-files from Delft3D-FM also called Dflow (operational version)" begin
      include("test_dflow.jl")
   end

   @testset "Read and interpolate map-files from Delft3D-FM (draft version including 3D flow)" begin
   include("test_delft3dfm.jl")
   end

   @testset "Convert delft3d_fm output to zarr" begin
   include("test_delft3dfm_to_zarr.jl")
   end

   @testset "Read and interpolate ERA5 files " begin
      include("test_era5.jl")
   end
      
   @testset "Read and interpolate Zarr files " begin
      include("test_zarr_grid.jl")
   end
      
   @testset "Read and interpolate NetCDF grid files from Matroos" begin
      include("test_matroos_grid.jl")
   end

   @testset "Read and interpolate NetCDF grid files from CMEMS" begin
      include("test_cmems_grid.jl")
   end

   @testset "Particle core" begin
      include("test_particles_core.jl")
   end
    
end

#clear cache
if isdir(joinpath(pwd(),".cache"))
   rm(joinpath(pwd(),".cache"),recursive=true)
end
