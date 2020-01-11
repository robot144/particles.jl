using Test, Particles, Printf, NetCDF

#clear cache
if isdir(joinpath(pwd(),".cache"))
   rm(joinpath(pwd(),".cache"),recursive=true)
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

@testset "Read and interpolate map-files from Delft3D-FM also called Dflow" begin
   include("test_dflow.jl")
end

@testset "Read and interpolate ERA5 files " begin
   include("test_era5.jl")
end

@testset "Read and interpolate NetCDF grid files from Matroos" begin
   include("test_matroos_grid.jl")
end

@testset "Particle core" begin
   include("test_particles_core.jl")
end

#clear cache
if isdir(joinpath(pwd(),".cache"))
   rm(joinpath(pwd(),".cache"),recursive=true)
end
