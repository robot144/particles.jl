#
# testing cartesian_grid.jl 
#
using NetCDF

#
# low-level tests
#

function test1() #CartesianXYGrid 
   xgrid=collect(range(0.0,1.0,step=0.1))
   ygrid=collect(range(0.0,2.0,step=0.1))
   grid=CartesianXYGrid(xgrid,ygrid)
   @test grid.bbox==[0.0,1.0,0.0,2.0]
   @test length(grid.xnodes)==11
   @test length(grid.ynodes)==21
   println("=== start: grid ===")
   dump(grid)
   println("===  end: grid ===")
   @test in_bbox(grid,0.5,0.5)==true
   @test in_bbox(grid,1.5,0.5)==false
   #spatial interpolation
   z=[ sin(x)*cos(y) for x in grid.xnodes , y in grid.ynodes ]
   x1=0.01
   y1=0.01
   z1_true=sin(x1)*cos(y1)
   xi1,yi1 = find_index(grid,x1,y1)
   @test (xi1,yi1)==(1,1)
   m1,n1,w1=find_index_and_weights(grid,x1,y1)
   @test m1[1]==1
   @test length(m1)==4
   @test length(n1)==4
   @test length(w1)==4
   @test w1[1]≈0.9*0.9 #bilinear
   z1=apply_index_and_weights(m1,n1,w1,z,NaN)
   @test abs(z1-z1_true)<1e-3
   zi=interpolate(grid,x1,y1,z,NaN)
   @test abs(zi-z1_true)<1e-3
   # point 2 is out of bounds
   x2=0.01
   y2=3.01
   xi2,yi2 = find_index(grid,x2,y2)
   @test (xi2,yi2)==(-1,-1)
   m2,n2,w2=find_index_and_weights(grid,x2,y2)
   @test m2[2]==-1
   z2=apply_index_and_weights(m2,n2,w2,z,NaN)
   @test isnan(z2)
end

function test2() #Cartesian space-time interpolation
   xgrid=collect(range(0.0,1.0,step=0.1))
   ygrid=collect(range(0.0,2.0,step=0.1))
   grid=CartesianXYGrid(xgrid,ygrid)
   times=collect(range(0.0,3.0,step=0.1))
   pressure(x,y,t)=1.0*x+10.0*y+100.0*t
   p=[pressure(x,y,t) for x in xgrid, y in ygrid, t in times]

   xyt=CartesianXYTGrid(grid,times,p,"p",0.0)
   @test xyt.name=="p"
   @test xyt.scaling==1.0
   @test xyt.offset==0.0
   @test xyt.missing_value==0.0
   @test length(xyt.cache)==3
   @test size(xyt.cache[3])==(11,21)
   @test all(xyt.time_cache.≈[0.0,0.1,0.2])
   @test xyt.time_cache_index==3

   #low level first
   update_cache(xyt,0.1) #no change to cache
   @test xyt.time_cache_index==3
   update_cache(xyt,0.3) #advance to next time
   @test xyt.time_cache_index==4
   update_cache(xyt,0.4) #advance to next time
   @test xyt.time_cache_index==5
   update_cache(xyt,1.0) #big step forward - refresh cache
   @test xyt.time_cache_index==13
   update_cache(xyt,1.51) #big step forward - refresh cache
   @test xyt.time_cache_index==18
   p1=get_map_slice(xyt,1)
   @test size(p1)==(11,21)
   @test p1[2,1]≈0.1
   @test p1[1,2]≈1.0
   update_cache(xyt,0.0) #no change to cache
   w=weights(xyt,0.025)
   @test all(w.≈(0.75,0.25,0.0))
   w=weights(xyt,0.125)
   @test all(w.≈(0.0,0.75,0.25))

   #higher level interpolation
   pi1=interpolate(xyt,0.1,0.0,0.0,999.0)
   @test pi1≈0.1
   pi2=interpolate(xyt,0.0,0.1,0.0,999.0)
   @test pi2≈1.0
   pi3=interpolate(xyt,0.1,0.0,0.1,999.0)
   @test pi3≈10.1
end

function test3() #CartesianXYGrid with one dimension of length 1
   xgrid=collect(range(0.0,1.0,step=0.1))
   ygrid=[0.0]
   grid=CartesianXYGrid(xgrid,ygrid)
   @test grid.bbox==[0.0,1.0,0.0,0.0]
   @test length(grid.xnodes)==11
   @test length(grid.ynodes)==1
   println("=== start: grid ===")
   dump(grid)
   println("===  end: grid ===")
   @test in_bbox(grid,0.5,0.0)==true
   @test in_bbox(grid,1.5,0.0)==false
   #spatial interpolation
   z=[ sin(x)*cos(y) for x in grid.xnodes , y in grid.ynodes ]
   x1=0.01
   y1=0.0
   z1_true=sin(x1)*cos(y1)
   xi1,yi1 = find_index(grid,x1,y1)
   @test (xi1,yi1)==(1,1)
   m1,n1,w1=find_index_and_weights(grid,x1,y1)
   println("m1=$(m1), n1=$(n1), w1=$(w1)")
   @test m1[1]==1
   @test length(m1)==4
   @test length(n1)==4
   @test length(w1)==4
   @test w1[1]≈0.9 #linear
   z1=apply_index_and_weights(m1,n1,w1,z,NaN)
   @test abs(z1-z1_true)<1e-3
   zi=interpolate(grid,x1,y1,z,NaN)
   @test abs(zi-z1_true)<1e-3
   # point 2 is out of bounds
   x2=1.01 #>1.0
   y2=0.0
   xi2,yi2 = find_index(grid,x2,y2)
   @test (xi2,yi2)==(-1,-1)
   m2,n2,w2=find_index_and_weights(grid,x2,y2)
   @test m2[2]==-1
   z2=apply_index_and_weights(m2,n2,w2,z,NaN)
   @test isnan(z2)
end

function test4() #Cartesian space-time interpolation with one dimension of length 1
   xgrid=collect(range(0.0,1.0,step=0.1))
   ygrid=[0.0]
   grid=CartesianXYGrid(xgrid,ygrid)
   times=collect(range(0.0,3.0,step=0.1))
   pressure(x,y,t)=1.0*x+10.0*y+100.0*t
   p=[pressure(x,y,t) for x in xgrid, y in ygrid, t in times]

   xyt=CartesianXYTGrid(grid,times,p,"p",0.0)
   @test xyt.name=="p"
   @test xyt.scaling==1.0
   @test xyt.offset==0.0
   @test xyt.missing_value==0.0
   @test length(xyt.cache)==3
   @test size(xyt.cache[3])==(11,1)
   @test all(xyt.time_cache.≈[0.0,0.1,0.2])
   @test xyt.time_cache_index==3

   #low level first
   update_cache(xyt,0.1) #no change to cache
   @test xyt.time_cache_index==3
   update_cache(xyt,0.3) #advance to next time
   @test xyt.time_cache_index==4
   update_cache(xyt,0.4) #advance to next time
   @test xyt.time_cache_index==5
   update_cache(xyt,1.0) #big step forward - refresh cache
   @test xyt.time_cache_index==13
   update_cache(xyt,1.51) #big step forward - refresh cache
   @test xyt.time_cache_index==18
   p1=get_map_slice(xyt,1)
   println("p1=$(p1)")
   @test size(p1)==(11,1)
   @test isnan(p1[1,1]) #0.0 is missing value
   @test p1[2,1]≈0.1
   @test p1[3,1]≈0.2

   update_cache(xyt,0.0) #no change to cache
   w=weights(xyt,0.025)
   @test all(w.≈(0.75,0.25,0.0))
   w=weights(xyt,0.125)
   @test all(w.≈(0.0,0.75,0.25))

   #higher level interpolation
   pi1=interpolate(xyt,0.1,0.0,0.0,999.0)
   @test pi1≈0.1
   pi2=interpolate(xyt,0.2,0.0,0.0,999.0)
   @test pi2≈0.2
   pi3=interpolate(xyt,0.1,0.0,0.1,999.0)
   @test pi3≈10.1
end

test1()
test2()
# similar to test2 but with one dimension of length 1
test3()
test4()


