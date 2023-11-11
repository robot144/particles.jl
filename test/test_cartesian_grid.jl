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

function test5() #CartesianXYZGrid with layers for z no times
   xgrid=collect(range(0.0,1.0,step=0.2))
   ygrid=collect(range(0.0,2.0,step=0.2))
   nlayers=3
   zgrid=zeros(length(xgrid),length(ygrid),nlayers+1) #depths of layer interfaces
   values=zeros(length(xgrid),length(ygrid),nlayers)  #values to be interpolated at each layer
   for k in 1:nlayers
      for j in 1:length(ygrid)
         for i in 1:length(xgrid)
            x = xgrid[i]
            y = ygrid[j]
            d = 10.0 + 5.0*exp(-((x-0.5)^2+(y-1.0)^2)/0.1) #depth
            zgrid[i,j,k+1]= -d*k/nlayers
            values[i,j,k]= -d*(k-0.5)/nlayers # z value at layer centre
         end
      end
   end
   grid=CartesianXYZGrid(xgrid,ygrid,zgrid)
   # xy_grid is now in grid.xy_grid
   @test grid.xy_grid.bbox==[0.0,1.0,0.0,2.0]
   @test length(grid.xy_grid.xnodes)==6
   @test length(grid.xy_grid.ynodes)==11
   println("=== start: grid ===")
   dump(grid)
   println("===  end: grid ===")
   @test in_bbox(grid,0.5,0.5)==true
   @test in_bbox(grid,1.5,0.5)==false
   #interpolate at a point
   x0=0.5;y0=0.5;z0=-2.0
   # double check in x,y space
   xi,yi=find_index(grid.xy_grid,x0,y0)
   @test (xi,yi)==(3,3)
   println("x[xi]=$(xgrid[xi]) y[yi]=$(ygrid[yi])")
   (x_indices,y_indices,weights)=find_index_and_weights(grid.xy_grid,x0,y0)
   println("x_indices=$(x_indices) y_indices=$(y_indices) weights=$(weights)")
   # now with z
   val_interp0 = interpolate(grid,x0,y0,z0,values,999.0)
   println("val_interp=$(val_interp0), z=$(z0)")
   @test val_interp0≈z0
   # point outside grid
   x1=1.5;y1=0.5;z1=-3.0
   val_interp1 = interpolate(grid,x1,y1,z1,values,999.0)
   println("val_interp=$(val_interp1), z=$(z1)")
   @test val_interp1≈999.0
   # point above surface
   x2=0.5;y2=0.5;z2=3.0
   val_interp2 = interpolate(grid,x2,y2,z2,values,999.0)
   println("val_interp=$(val_interp2), z=$(z2)")
   @test val_interp2≈999.0
end

function test6() #Cartesian space-time interpolation with 3 dimensions plus time
   xgrid=collect(range(0.0,1.0,step=0.1))
   ygrid=collect(range(0.0,2.0,step=0.2))
   times=collect(range(0.0,3.0,step=0.1))
   nlayers=3
   depth(x,y)=10.0 + 5.0*exp(-((x-0.5)^2+(y-1.0)^2)/0.1) #depth
   waterlevel(x,y,t)=0.01*t
   layer_ifaces(x,y,k,t)= waterlevel(x,y,t)-(waterlevel(x,y,t)-depth(x,y))*k/nlayers
   zgrid=[layer_ifaces(x,y,k,t) for x in xgrid, y in ygrid, k in 0:nlayers, t in times]
   pressure(x,y,k,t)= -(waterlevel(x,y,t)-depth(x,y))*(k-0.5)/nlayers # z value at layer centre
   p=[pressure(x,y,k,t) for x in xgrid, y in ygrid, k in 1:nlayers, t in times]

   grid=CartesianXYZGrid(xgrid,ygrid,[]) #zgrid will come later
   xyzt=CartesianXYZTGrid(grid,times,p,zgrid,"p",0.0,9999.0)
   @test xyzt.name=="p"
   @test xyzt.scaling==1.0
   @test xyzt.offset==0.0
   @test xyzt.missing_value==0.0
   @test length(xyzt.cache)==3
   @test size(xyzt.cache[3])==(11,11,3)
   @test all(xyzt.time_cache.≈[0.0,0.1,0.2])
   @test xyzt.time_cache_index==3

   #low level first
   update_cache(xyzt,0.1) #no change to cache
   @test xyzt.time_cache_index==3
   update_cache(xyzt,0.3) #advance to next time
   @test xyzt.time_cache_index==4
   update_cache(xyzt,0.4) #advance to next time
   @test xyzt.time_cache_index==5
   update_cache(xyzt,1.0) #big step forward - refresh cache
   @test xyzt.time_cache_index==13
   update_cache(xyzt,1.51) #big step forward - refresh cache
   @test xyzt.time_cache_index==18
   p1=get_map_slice(xyzt,1)
   println("p1=$(p1)")
   @test size(p1)==(11,11,3)
   @test p1[1,1,1]≈1.6666697722109767
   @test p1[2,1,1]≈1.6666743050731136
   @test p1[1,3,1]≈1.6685357230995714
   @test p1[1,1,3]≈8.333348861054883
   z1=get_zmap_slice(xyzt,1)
   println("z1=$(z1)")

   update_cache(xyzt,0.0) #no change to cache
   w=weights(xyzt,0.025)
   @test all(w.≈(0.75,0.25,0.0))
   w=weights(xyzt,0.125)
   @test all(w.≈(0.0,0.75,0.25))

   #higher level interpolation
   x1=0.1;y1=0.0;z1=0.0;t1=0.0
   pi1=interpolate(xyzt,x1,y1,z1,t1,999.0)
   @test pi1≈999.0 #pressure(x1,y1,1,t1)
   x2=0.2;y2=0.2;z2=2.0;t2=0.0
   pi2=interpolate(xyzt,x2,y2,z2,t2,999.0)
   @test pi2≈2.0
   x3=0.4;y3=0.6;z3=8.0;t3=0.0
   pi3=interpolate(xyzt,x3,y3,z3,t3,999.0)
   @test pi3≈8.0
   x4=0.4;y4=0.6;z4=8.0;t4=0.0
   pi4=interpolate(xyzt,x4,y4,z4,t4,999.0)
   @test pi4≈8.0
end

function test7() #Cartesian space-time interpolation with 3 dimensions plus time - values at layer interfaces
   xgrid=collect(range(0.0,1.0,step=0.1))
   ygrid=collect(range(0.0,2.0,step=0.2))
   times=collect(range(0.0,3.0,step=0.1))
   nlayers=3
   depth(x,y)=10.0 + 5.0*exp(-((x-0.5)^2+(y-1.0)^2)/0.1) #depth
   waterlevel(x,y,t)=0.01*t
   layer_ifaces(x,y,k,t)= waterlevel(x,y,t)-(waterlevel(x,y,t)-depth(x,y))*k/nlayers
   zgrid=[layer_ifaces(x,y,k,t) for x in xgrid, y in ygrid, k in 0:nlayers, t in times]
   pressure(x,y,k,t)= -(waterlevel(x,y,t)-depth(x,y))*k/nlayers # z value at layer interface
   p=[pressure(x,y,k,t) for x in xgrid, y in ygrid, k in 0:nlayers, t in times]

   spherical=false
   value_at_iface=true
   grid=CartesianXYZGrid(xgrid,ygrid,[],spherical,value_at_iface) #zgrid will come later
   xyzt=CartesianXYZTGrid(grid,times,p,zgrid,"p",0.0,9999.0)
   @test xyzt.name=="p"
   @test xyzt.scaling==1.0
   @test xyzt.offset==0.0
   @test xyzt.missing_value==0.0
   @test length(xyzt.cache)==3
   @test size(xyzt.cache[3])==(11,11,4)
   @test all(xyzt.time_cache.≈[0.0,0.1,0.2])
   @test xyzt.time_cache_index==3

   #higher level interpolation
   x1=0.1;y1=0.0;z1=0.0;t1=0.0
   pi1=interpolate(xyzt,x1,y1,z1,t1,999.0)
   @test pi1≈999.0 #pressure(x1,y1,1,t1)
   x2=0.2;y2=0.2;z2=2.0;t2=0.0
   pi2=interpolate(xyzt,x2,y2,z2,t2,999.0)
   @test pi2≈999.0
   x3=0.4;y3=0.6;z3=8.0;t3=0.0
   pi3=interpolate(xyzt,x3,y3,z3,t3,999.0)
   @test pi3≈8.0
   x4=0.4;y4=0.6;z4=8.0;t4=0.0
   pi4=interpolate(xyzt,x4,y4,z4,t4,999.0)
   @test pi4≈8.0
end

test1()
test2()
# similar to test2 but with one dimension of length 1
test3()
test4()
test5()
test6()
test7()

