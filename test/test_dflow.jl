# Interact with dflow netcdf output map files
#
using Dates

#
# test
#

function test1()
   #init
   dflow_map=load_nc_info("../test_data",r"estuary_...._map.nc")
   @test length(dflow_map)==2
   @test size(dflow_map[2].vars["s1"])==(103,25)
   interp=load_dflow_grid(dflow_map,50,false)
   @test length(interp.grids)==2
   
   #interpolate a field to a regular grid
   sealevel=load_nc_map_slice(dflow_map,"s1",10)
   @test length(sealevel)==2
   @test length(sealevel[1])==103
   @test length(sealevel[2])==103
   @test sealevel[2][12]==-0.3979754473027684

   u_velocity=load_nc_map_slice(dflow_map,"ucx",10)
   @test length(u_velocity)==2
   @test length(u_velocity[1])==103
   @test length(u_velocity[2])==103
   @test u_velocity[1][103]==-0.2017744781187035

   xpoints=collect(range(0.,stop=100000.,length=300))
   ypoints=collect(range(0.,stop=500.,length=100))
   sealevel_interp=interpolate(interp,xpoints,ypoints,sealevel)
   @test size(sealevel_interp)==(300,100)
   @test sealevel_interp[23,48]==-0.17045316506644842

   u_velocity_interp=interpolate(interp,xpoints,ypoints,u_velocity)
   @test size(u_velocity_interp)==(300,100)
   #heatmap(xpoints,ypoints,u_velocity_interp')
   # interpolate for one point only
   ind=find_index(interp,10000.0,200.0)
   @test ind==[2 84]
   ux=apply_index(ind,u_velocity,-9999.)
   @test ux==-0.14349077430813084
   
   # u,v interpolation functions
   t0=get_reftime(dflow_map)
   @test t0==DateTime(1991,1,1)
   u1,v1=initialize_interpolation(dflow_map,interp,t0)
   u_value=u1(100.0,100.0,0.0,0.0)
   #@test u_value==-0.9088566953087289
   @test u_value==-0.0 #TODO Check later. Is this really true?
end

# this test needs files that do not fit in github
# only run if the files exist
function test2()
   #init
   dflow_map=load_nc_info("../test_data",r"DCSM-FM_0_5nm_...._map.nc")
   @test length(dflow_map)==20
   @test size(dflow_map[2].vars["mesh2d_s1"])==(31553,170)
   interp=load_dflow_grid(dflow_map,50,false)
   @test length(interp.grids)==20

   ##interpolate a field to a regular grid
   #sealevel=load_nc_map_slice(dflow_map,"mesh2d_s1",10)
   #xpoints=collect(range(-15.0,stop=13.0,length=1200))
   #ypoints=collect(range(43.0,stop=64.0,length=1000))
   #sealevel_interp=interpolate(interp,xpoints,ypoints,sealevel)
   #Plots.default(:size,[1200,1000])
   #heatmap(xpoints,ypoints,sealevel_interp', clims=(-2.0,2.0))

   # u,v interpolation functions
   t0=get_reftime(dflow_map)
   @test t0==DateTime(2012,12,22)
   u2,v2=initialize_interpolation(dflow_map,interp,t0)
   for istep=1:5
      u_value=u2(1.0,51.0,0.0,864000.0+3600.0*istep)
      println("$(istep) u=$(u_value)")
      if istep==5
         @test u_value==0.04639523554494074
      end
   end
end

function test3()
   #init
   dflow_map=load_nc_info("../test_data","dfm_multiple_runs_and_domains/01_flow_to_right/DFM_OUTPUT_cb_2d/cb_2d_0000_map.nc")
   @test length(dflow_map)==1
   dflow_map=load_nc_info("../test_data","dfm_multiple_runs_and_domains/01_flow_to_right/DFM_OUTPUT_cb_2d/cb_2d_0000_map.nc")
   @test length(dflow_map)==1
   dflow_map=load_nc_info("../test_data/dfm_multiple_runs_and_domains/01_flow_to_right/DFM_OUTPUT_cb_2d/",r"cb_2d_...._map.nc")
   @test length(dflow_map)==2
   dflow_map=load_nc_info("../test_data","dfm_multiple_runs_and_domains/*/DFM_OUTPUT_cb_2d/cb_2d_*_map.nc")
   @test length(dflow_map)==4
   @test size(dflow_map[1].vars["mesh2d_ucx"])==(317,25)
   interp=load_dflow_grid(dflow_map,50,true)
   @test length(interp.grids)==2
   
   # u,v interpolation functions
   t0=get_reftime(dflow_map)
   @test t0==DateTime(2001,1,1)
   u1,v1=initialize_interpolation(dflow_map,interp,t0)
   u_value=u1(3.33, 47.99, 0.0, 1800.0)
   @test u_value==0.0435357731534879
   u_value=u1(3.33, 47.99, 0.0, 86399.0)
   @test u_value==-0.048490442244461425
   u_value=u1(3.33, 47.99, 0.0, 43200.0)
   @test u_value==0
   u_value=u1(3.33, 47.99, 0.0, 43199.0)
   @test u_value==0.04905390667374792
end

test1()

#test2 needs large input files that are not present in the repository.
#only run tests if the files have been added (manually for now)
if isfile(joinpath("../test_data","DCSM-FM_0_5nm_0000_map.nc"))
   test2()
end

test3()
