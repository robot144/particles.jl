#
# testing unstructured_grid.jl 
#
using NetCDF

#
# low-level tests
#

function test1()
   xpol1=[1.,0.,-1]
   ypol1=[-1.,1.,-1.]
   wn1_in=winding_number(0.,0.,xpol1,ypol1)
   @test wn1_in==1
   wn1_out=winding_number(1.,0.,xpol1,ypol1)
   @test wn1_out==0
   
   xnodes=[0.0,1.0,2.0,0.0,1.0,2.0,0.0,1.0,2.0,0.0]
   ynodes=[0.0,0.0,0.0,1.0,1.0,1.0,2.0,2.0,2.0,3.0]
   edges=[1 2 5 4;2 3 6 5;4 5 8 7;5 6 9 8;7 8 10 -1]
   nnodes=nodes_per_cell(edges)
   @test length(nnodes)==5
   @test nnodes[1]==4
   @test nnodes[4]==4
   @test nnodes[5]==3
   icell=find_first_cell(0.5,0.5,xnodes,ynodes,edges,nnodes)
   @test icell==1
   icell=find_first_cell(1.5,0.5,xnodes,ynodes,edges,nnodes)
   @test icell==2
   icell=find_first_cell(0.25,2.5,xnodes,ynodes,edges,nnodes)
   @test icell==5
   icell=find_first_cell(-0.5,0.5,xnodes,ynodes,edges,nnodes)
   @test icell==-1
end

#
# medium-level tests
#

function test2()
   xnodes=[0.0,1.0,2.0,0.0,1.0,2.0,0.0,1.0,2.0,0.0,3.0]
   ynodes=[0.0,0.0,0.0,1.0,1.0,1.0,2.0,2.0,2.0,3.0,0.0]
   edges=[1 2 5 4;2 3 6 5;4 5 8 7;5 6 9 8;7 8 10 -1; 3 11 6 -1]
   #grid objects
   mygrid=Grid(xnodes,ynodes,edges,0) #turn off indexing
   icell=find_first_cell(0.5,0.5,mygrid)
   @test icell==1
   icell=find_first_cell(1.5,0.5,mygrid)
   @test icell==2
   icell=find_first_cell(0.25,2.5,mygrid)
   @test icell==5
   icell=find_first_cell(-0.5,0.5,mygrid)
   @test icell==-1
   
   xpoints=collect(-1.0:0.01:3.0)
   ypoints=collect(-1.0:0.01:4.0)
   cells=find_cells(xpoints,ypoints,mygrid)
   (m,n)=size(cells)
   @test m==401
   @test n==501
   @test cells[200,200]==1

   dump(mygrid)
   println("=== start: create_node_tree ===")
   create_node_tree!(mygrid,2)
   println("===  end: create_node_tree ===")
   dump(mygrid)

   #now test with indexing
   icell=find_first_cell(0.5,0.5,mygrid)
   @test icell==1
   icell=find_first_cell(1.5,0.5,mygrid)
   @test icell==2
   icell=find_first_cell(0.25,2.5,mygrid)
   @test icell==5
   icell=find_first_cell(-0.5,0.5,mygrid)
   @test icell==-1

   cells=find_cells(xpoints,ypoints,mygrid)
   (m,n)=size(cells)
   @test m==401
   @test n==501
   @printf("cells[200,200]=%d\n",cells[200,200])
   @test cells[200,200]==1
end

function test3()
   #medium level interpolation on data from netcdf
   map0=ncinfo("../test_data/estuary_0000_map.nc")
   edges0=map0.vars["NetElemNode"][:,:]'
   xnodes0=map0.vars["NetNode_x"][:]
   ynodes0=map0.vars["NetNode_y"][:]
   grid0=Grid(xnodes0,ynodes0,edges0,10)
   dump(grid0)
   (m,n)=size(grid0.edges)
   @test m==103
   @test n==4
   print("grid0.bbox ")
   println(grid0.bbox)
   @test abs(grid0.bbox[1]-47501.0)<1e-3
   @test abs(grid0.bbox[2]-99001.0)<1e-3
   @test abs(grid0.bbox[3]-1.)<1e-3
   @test abs(grid0.bbox[4]-501.)<1e-3
   
   itime=25 #last timestep in this example
   waterlevel0=map0.vars["s1"][:,itime]
   xpoints0=range(grid0.bbox[1],stop=grid0.bbox[2],length=100)
   ypoints0=range(grid0.bbox[3],stop=grid0.bbox[4],length=10)
   xpoint=50000.
   ypoint=250.
   cell0=find_first_cell(xpoint,ypoint,grid0)
   @test cell0==99
   cells0=find_cells(xpoints0,ypoints0,grid0)
   wl_interp0=zeros(Float64,size(cells0))
   get_values_by_cells!(wl_interp0,cells0,waterlevel0)
   @printf("maximum(wl_interp0)=%f\n",maximum(wl_interp0))
   @test abs(maximum(wl_interp0)-0.12916)<1e-3
   #return (xpoints0,ypoints0,wl_interp0)
end

#
# high-level tests
#

function test4(nmin=50)
   map=[]
   push!(map,ncinfo("../test_data/estuary_0000_map.nc"))
   push!(map,ncinfo("../test_data/estuary_0001_map.nc"))
   interp=Interpolator()
   for i=1:length(map)
      edges_temp=map[i].vars["NetElemNode"][:,:]'
      xnodes_temp=map[i].vars["NetNode_x"][:]
      ynodes_temp=map[i].vars["NetNode_y"][:]
      @printf("- index computation\n")
      @time grid_temp=Grid(xnodes_temp,ynodes_temp,edges_temp,nmin)
      dump(grid_temp)
      add_grid!(interp,grid_temp)
   end
   xpoints=range(0.,stop=100000,length=2000)
   ypoints=range(0,stop=500,length=2000)
   xpoints2=range(0.,stop=100000,length=2001) #create new points that are not already cached
   ypoints2=range(0,stop=500,step=2001)
   itime=25 #last timestep in this example
   waterlevel0=map[1].vars["s1"][:,itime]
   waterlevel1=map[2].vars["s1"][:,itime]
   waterlevel=[]
   push!(waterlevel,waterlevel0)
   push!(waterlevel,waterlevel1)
   @printf("1 - compute new weights\n")
   @time wl_interp=interpolate(interp,xpoints,ypoints,waterlevel)
   @printf("2 - use existing weights\n")
   @time wl_interp=interpolate(interp,xpoints,ypoints,waterlevel)
   #@printf("3 - compute new weights\n")
   #@time wl_interp2=interpolate(interp,xpoints2,ypoints2,waterlevel)
   #@printf("4 - use existing weights\n")
   #@time wl_interp2=interpolate(interp,xpoints2,ypoints2,waterlevel)
   (m,n)=size(wl_interp)
   @test m==length(xpoints)
   @test n==length(ypoints)
   #@test abs(maximum(wl_interp)-0.6241)<1e-3
   #@test abs(minimum(wl_interp)+0.3729)<1e-3
   #pcolormesh(xpoints,ypoints,wl_interp')
end

function test5(nmin)
   #large test for measuring performance
   map=[]
   push!(map,ncinfo("../output/gtsm_fine_0000_map.nc"))
   push!(map,ncinfo("../output/gtsm_fine_0001_map.nc"))
   push!(map,ncinfo("../output/gtsm_fine_0002_map.nc"))
   push!(map,ncinfo("../output/gtsm_fine_0003_map.nc"))
   interp=Interpolator()
   for i=1:length(map)
      edges_temp=map[i].vars["NetElemNode"][:,:]'
      xnodes_temp=map[i].vars["NetNode_x"][:]
      ynodes_temp=map[i].vars["NetNode_y"][:]
      @printf("- index computation\n")
      @time grid_temp=Grid(xnodes_temp,ynodes_temp,edges_temp,nmin,true)
      #dump(grid_temp)
      add_grid!(interp,grid_temp)
   end
   xpoints=collect(range(-180.,stop=180.,length=300))
   ypoints=collect(range(-90.,stop=90.,length=300))
   #xpoints=range(-18.,stop=18.,length=1000)
   #ypoints=range(40.,stop=70.,length=1000)
   itime=25 #last timestep in this example
   waterlevel=[]
   for i=1:length(map)
      push!(waterlevel,map[i].vars["s1"][:,itime])
   end
   println("interpolation")
   @time wl_interp=interpolate(interp,xpoints,ypoints,waterlevel)
   @time wl_interp=interpolate(interp,xpoints,ypoints,waterlevel)
   return (xpoints,ypoints,wl_interp)
end

test1()
test2()
test3()
test4(20)

#This is just for debugging and performance tests
#(x5,y5,z5)=test5(10)
#using Plots
#contourf(5,y5,z5',clims=(-1.0,1.0))
#Profile.clear_malloc_data()
#@time test5(10)
#pcolormesh(xpoints,ypoints,wl_interp',vmin=-1,vmax=1)
