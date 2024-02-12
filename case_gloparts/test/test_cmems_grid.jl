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

function test2()
   #init
   cmems_data = CmemsData("../test_data","u0_2021-10-29_00-00-00_2021-11-07_00-00-00.nc")
   
   u1=load_map_slice(cmems_data,"uo",1) # u-velocity
   @test size(u1)==(733,277)
   @test isnan(u1[1,1])
   @test isapprox(u1[100,100],0.21790215,atol=0.001)

   t0=get_reftime(cmems_data)
   @test t0==DateTime(2021,10,29)

   times=get_times(cmems_data,t0)
   @test times[1]≈1800.0
   @test times[2]≈5400.0

   t2=as_DateTime(cmems_data,t0,times[2])
   @test t2==DateTime("2021-10-29T01:30:00")


   u=initialize_interpolation(cmems_data,"uo",t0,NaN) #water velocity x-dir
   u1=u(-73.0,35.0,0.0,1800.0)
   @test u1≈0.04947685987763097
   u2=u(-73.0,35.0,0.0,3600.0)
   @test u2≈0.037575478067336805
end

function test3()
   #init
   cmems_datas = CmemsData("../test_data", r"cmems_map_u_.+_time..nc"; lon = "x", lat = "y") # Regex
   @test length(cmems_datas) == 2
   cmems_datas = CmemsData("../test_data", "cmems_map_u_*_time*.nc"; lon = "x", lat = "y") # Glob
   @test length(cmems_datas) == 2
   cmems_data = cmems_datas[1]

   t0 = get_reftime(cmems_data)
   @test t0 == DateTime(2022,05,20)

   times = get_times(cmems_data,t0)
   @test times[1] ≈ 0.0
   @test times[2] ≈ 3600.0

   t2 = as_DateTime(cmems_data, t0, times[2])
   @test t2 == DateTime("2022-05-20T01:00:00")


   u = initialize_interpolation(cmems_datas, "uo", t0, NaN) #water velocity x-dir
   u1 = u(0.5, 52, 0.0, 1800.0) # only available in 'time1'-file
   @test u1 ≈ -1
   u1 = u(0.5, 52, 0.0, 21600.0) # available in both 'time1'- and 'time2'-file (so pick 'time2'-file)
   @test u1 ≈ 1
   u1 = u(0.5, 52, 0.0, 90000.0) # only available in 'time2'-file
   @test u1 ≈ 1
end

test1()
test2()
test3()
