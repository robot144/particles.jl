# tests for main Particles core routines from particles_core.jl

exampledir=joinpath(pwd(),"examples")
if !isdir(exampledir)
   exampledir=joinpath("..","examples")
end
testdatadir=joinpath(pwd(),"test_data")
if !isdir(testdatadir)
   testdatadir=joinpath("..","test_data")
end

#test 1
include(joinpath(exampledir,"example_loop.jl"))
@test d["tstart"]==0.0
@test d["tend"]==2.5
d["plot_maps"]=false #turn off plotting in unit tests
@eval run_simulation(d)
#look at output    TODO: still taking one timestep too much due to rounding
@test d["particles"][4,1]≈2.50

#test 2
include(joinpath(exampledir,"example_2dv_sideview.jl"))
@test d["tstart"]==0.0
@test d["tend"]==90.0
d["plot_maps"]=false #turn off plotting in unit tests
@eval run_simulation(d)
@test d["particles"][1,1]≈90.0

#test 3
include(joinpath(exampledir,"example_dflow_2d_estuary.jl"))
@test d["tstart"]==0.0
@test d["tend"]≈(24*3600.0)
d["plot_maps"]=false #turn off plotting in unit tests
@eval run_simulation(d)
@test d["particles"][4,1]≈(24*3600.0)

#test 4
if isfile(joinpath(testdatadir,"DCSM-FM_0_5nm_0019_map.nc"))
   include(joinpath(exampledir,"example_dflow_2d_dcsm.jl"))
   @test d["tstart"]≈(240*3600.0)
   @test d["tend"]≈(408*3600.0)
   d["plot_maps"]=false #turn off plotting in unit tests
   @eval run_simulation(d)
   @test d["particles"][3,1]≈(168*3600.0)
else
   println("No test data found for test example_dflow_2d_dcsm.jl")
end

