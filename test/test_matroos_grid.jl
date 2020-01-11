# Interact with dflow netcdf output map files
#
using Dates

#
# test
#

function test1()
   #init
   matroos_data = MatroosData("../test_data","swan_dcsm_20191025.nc")
   
   hs=load_map_slice(matroos_data,"wave_height_hm0",1)
   @test size(hs)==(421,481)
   @test hs[1,1]≈4.051f0

   t0=get_reftime(matroos_data)
   @test t0==DateTime(2019,10,25)

   times=get_times(matroos_data,t0)
   @test times[1]≈21600.0
   @test times[2]≈25200.0

   t2=as_DateTime(matroos_data,t0,times[2])
   @test t2==DateTime("2019-10-25T07:00:00")

   hs=initialize_interpolation(matroos_data,"wave_height_hm0",t0)
   hs1=hs(-11.975,47.98493499999999,0.0,21600.0)
   @test isapprox(hs1,4.022,atol=0.001)
   hs2=hs(-12.025,48.018265,0.0,21600.0)
   @test isapprox(hs2,4.097,atol=0.001)
end


test1()
