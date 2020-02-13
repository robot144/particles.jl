# Interact with CMEMS netcdf output map files
#
using Dates

#
# test
#

function test1()
   #init
   cmems_data = CmemsData("../test_data","cmems_20140301_06.nc")
   
   ice1=load_map_slice(cmems_data,"siconc",1) # sea-ice concentration
   @test size(ice1)==(481,181)
   @test isnan(ice1[1,1])
   @test isapprox(ice1[1,end],0.9666,atol=0.001)

   t0=get_reftime(cmems_data)
   @test t0==DateTime(2014,3,1)

   times=get_times(cmems_data,t0)
   @test times[1]≈43200.0
   @test times[2]≈129600.0

   t2=as_DateTime(cmems_data,t0,times[2])
   @test t2==DateTime("2014-03-02T12:00:00")


   u=initialize_interpolation(cmems_data,"uo",t0,NaN) #water velocity x-dir
   u1=u(20.0,75.0,0.0,43200.0)
   @test u1≈-0.04763681870922009
   u2=u(20.0,75.0,0.0,129600.0)
   @test u2≈-0.037204620036251405
   v=initialize_interpolation(cmems_data,"vo",t0,NaN) #water velocity y-dir
   v1=v(20.0,75.0,0.0,43200.0)
   @test v1≈0.04879613983597007
end


test1()
