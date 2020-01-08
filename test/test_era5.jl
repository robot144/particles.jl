# Interact with dflow netcdf output map files
#
using Dates

#
# test
#

function test1()
   #init
   era_data = EraData("../test_data","era5_wind_20140301_06.nc")
   
   p1=load_map_slice(era_data,"msl",1)
   @test size(p1)==(161,61)
   @test p1[1,1]≈102574.4477038286

   t0=get_reftime(era_data)
   @test t0==DateTime(2014,3,1)

   times=get_times(era_data,t0)
   @test times[1]≈0.0
   @test times[2]≈3600.0

   t2=as_DateTime(era_data,t0,times[2])
   @test t2==DateTime("2014-03-01T01:00:00")

   # u,v interpolation functions
   t0=get_reftime(era_data)
   @test t0==DateTime(2014,3,1)

   p=initialize_interpolation(era_data,"msl",t0)
   p1=p(20.0,75.0,0.0,0.0)
   @test p1≈101877.70372941876
   p2=p(20.0,75.0,0.0,3600.0)
   @test p2≈101926.02316008728
end


test1()
