# some tools for unstructured grids
using Printf
using Test
using NetCDF

mutable struct Grid
   #global list (by ref) even for subgrids
   xnodes::Array{Float64,1}     #x-coordinates of nodes
   ynodes::Array{Float64,1}     #y-coordinates of nodes
   edges::Array{Int64,2}        #list of node numbers around a cell. Row nr of entry is cell nr
                                #negative numbers for padding, eg [1,2,3,-1] for triangle 1,2,3
   #derived data for optimization
   edge_indices::Array{Int64,1} #in case of sugrids, the indices of edges at this level
   bbox::Array{Float64,1}       #bounding box [xmin xmax ymin ymax]
   nnodes::Array{Int64,1}       #cached count of nodes foreach cell

   subgrids::Array{Grid,1}      #multilevel storage for opimization
   function Grid(xnodes::Array,ynodes::Array,edges::AbstractMatrix,nmin=10,spherical::Bool=false)
      (nedges,maxNodesPerEdge)=size(edges)
      if(maxNodesPerEdge>6)
         if(nedges<=6)
            edges=collect(edges')
         else
            println("size(edges)=($(nedges),$(MaxNodesPerEdge))") 
         end
      end
      (nedges,maxNodesPerEdge)=size(edges)
      nnodes=nodes_per_cell(edges)
      bbox_temp=zeros(Float64,4)
      if spherical
         edge_indices=zeros(Int64,0)
         for i=1:nedges
	    compute_boundingbox!(bbox_temp,xnodes,ynodes,edges,nnodes,i)
            if (bbox_temp[2]-bbox_temp[1])<180.0
	       push!(edge_indices,i)
	    end
	 end
      else
         edge_indices=collect(1:nedges)
      end
      bbox=compute_boundingbox(xnodes,ynodes,edges,nnodes)
      this=new(xnodes,ynodes,edges,edge_indices,bbox,nnodes,[])
      if nmin>0
        create_node_tree!(this,nmin) 
      end
      return this
   end
   function Grid(xnodes::Array,ynodes::Array,edges::Array,edge_indices::Array,bbox::Array,nnodes::Array)
      new(xnodes,ynodes,edges,edge_indices,bbox,nnodes,[])
   end
end

mutable struct Interpolator
   grids::Array{Grid,1}  #grids for each domain
   outx_hash::UInt64    #hash code of previous output grid
   outy_hash::UInt64    #for caching the interpolation indices and weights
   cached_indices::Array{Any,1} #indices of previous interpolation
   function Interpolator()
      new([],hash([]),hash([]),[])
   end
end

#
# High-level functions
#

function add_grid!(interp::Interpolator,grid::Grid)
    push!(interp.grids,grid)
end

function get_bbox(interp::Interpolator)
   bbox=copy(interp.grids[1].bbox)
   for counter=2:length(interp.grids)
      temp_bbox=interp.grids[counter].bbox
      bbox[1]=min(bbox[1],temp_bbox[1]) #xmin
      bbox[2]=max(bbox[2],temp_bbox[2]) #xmax
      bbox[3]=min(bbox[3],temp_bbox[3]) #ymin
      bbox[4]=max(bbox[4],temp_bbox[4]) #ymax
   end
   return bbox
end

function interpolate(interp::Interpolator,xpoints,ypoints,invalues::Array)
   #check for mutiple or single domain
   if(typeof(invalues).parameters[1] <: Number)
      #input as single domain
      if length(interp.xnodes)>1
         error("Input values were given for single domain to multi-domain interpolator")
      end
      values=[invalues]
   end
   #compute hash and see if we can reuse interpolation weights
   outx_hash=hash(xpoints)
   outy_hash=hash(ypoints)
   if (interp.outx_hash!=outx_hash)||(interp.outy_hash!=outy_hash)
      #new output grid, so compute weights 
      println("compute new interpolation weights")
      interp.outx_hash=outx_hash
      interp.outy_hash=outy_hash
      cached_indices=[]
      for i=1:length(interp.grids)
         @printf("grid %d\n",i)
         @time cells=find_cells(xpoints,ypoints,interp.grids[i])
         push!(cached_indices,cells)
      end
      interp.cached_indices=cached_indices #store to cache for next call to this function
   else
      println("use cached interpolation weights")
      cached_indices=interp.cached_indices #retrieve from cache
   end
   outvalues=zeros(Float64,size(cached_indices[1]))
   for i=1:length(cached_indices) #loop over grids
      get_values_by_cells!(outvalues,cached_indices[i],invalues[i])
   end
   return outvalues
end

#
# low-level functions
# 

function nodes_per_cell(edges::AbstractMatrix)
   #each row of edges contains the node numbers of the cell 
   #padded with negative numbers
   #returns a vector with the number of nodes for each cell
   n,m = size(edges)
   nnodes=m*ones(Int64,n)
   for i=1:n
      for j=m:-1:3 #assume each cell has at least two nodes
         if(edges[i,j]<0)
            nnodes[i]=(j-1)
         end
      end
   end
   return nnodes
end

# compute winding number
# for simple counter clockwise polygons this gives 1 if point is inside and 0 outside
function winding_number(xpoint::AbstractFloat,ypoint::AbstractFloat,xspolygon::Array,yspolygon::Array)
   wn=0
   xprev=xspolygon[end]
   yprev=yspolygon[end]
   for i=1:length(xspolygon)
      xthis=xspolygon[i]
      ythis=yspolygon[i]
      # look for crossing of y=ypoint and x>xpoint, ie east of point
      if((yprev<ypoint)&&(ythis>=ypoint))
         #check for northward crossing east of point
	 s=(ypoint-yprev)/(ythis-yprev)
	 x=xprev+s*(xthis-xprev)
	 if x>xpoint
	    wn+=1
         end
      end
      if((yprev>=ypoint)&&(ythis<ypoint))
         #check for southward crossing east of point
	 s=(ypoint-yprev)/(ythis-yprev)
	 x=xprev+s*(xthis-xprev)
	 if x>xpoint
	    wn-=1
         end
      end
      xprev=xthis
      yprev=ythis
   end
   return wn
end

function find_first_cell(xpoint,ypoint,xnodes,ynodes,edges,nnodes)
   #for point (xpoint,ypoint) find cell with point insides
   #stops at first hit and return row number of edges array
   #return -1 if no cell is found
   #TODO start with linear search. Make more efficient later
   ##@printf("x=%f y=%f \n",xpoint,ypoint)
   m,n=size(edges)
   for i=1:m
      wn=0
      xprev=xnodes[edges[i,nnodes[i]]]
      yprev=ynodes[edges[i,nnodes[i]]]
      for j=1:nnodes[i]
         ##@printf("edge %d\n",j)
         xthis=xnodes[edges[i,j]]
         ythis=ynodes[edges[i,j]]
         ##@printf("(%f,%f) -> (%f,%f)\n",xprev,yprev,xthis,ythis)
         # look for crossing of y=ypoint and x>xpoint, ie east of point
         if((yprev<ypoint)&&(ythis>=ypoint))
            #check for northward crossing east of point
            s=(ypoint-yprev)/(ythis-yprev)
   	    x=xprev+s*(xthis-xprev)
   	    if x>xpoint
   	       wn+=1
               ##@printf("north\n")
            end
         end
         if((yprev>=ypoint)&&(ythis<ypoint))
            #check for southward crossing east of point
   	    s=(ypoint-yprev)/(ythis-yprev)
   	    x=xprev+s*(xthis-xprev)
   	    if x>xpoint
   	       wn-=1
               ##@printf("south\n")
            end
         end
         xprev=xthis
         yprev=ythis
      end
      #wn now contains winding number
      ##@printf("cell %d winding-number %d \n",i,wn)
      if wn!=0 
         return i
      end
   end
   return -1
end

function get_values_by_cells!(outvalues,cells,invalues)
   for i=1:length(cells)
      if(cells[i]>0)
         outvalues[i]=invalues[cells[i]]
      end
   end
end

function compute_boundingbox(xnodes,ynodes,edges,nnodes)
   #returns array with min(x),max(x),min(y),max(y)
   #only considers nodes that are part of a listed cell 
   xmin=xnodes[edges[1,1]]
   xmax=xnodes[edges[1,1]]
   ymin=ynodes[edges[1,1]]
   ymax=ynodes[edges[1,1]]
   m,n = size(edges)
   for i=1:m
   for j=1:nnodes[i]
      xmin=min(xmin,xnodes[edges[i,j]])
      xmax=max(xmax,xnodes[edges[i,j]])
      ymin=min(ymin,ynodes[edges[i,j]])
      ymax=max(ymax,ynodes[edges[i,j]])
   end
   end
   return [xmin,xmax,ymin,ymax]
end

function compute_boundingbox!(bbox,xnodes,ynodes,edges,nnodes,edge_index)
   #returns array with min(x),max(x),min(y),max(y)
   i=edge_index
   #only considers nodes that are part of a listed cell 
   xmin=xnodes[edges[i,1]]
   xmax=xnodes[edges[i,1]]
   ymin=ynodes[edges[i,1]]
   ymax=ynodes[edges[i,1]]
   for j=1:nnodes[i]
      xmin=min(xmin,xnodes[edges[i,j]])
      xmax=max(xmax,xnodes[edges[i,j]])
      ymin=min(ymin,ynodes[edges[i,j]])
      ymax=max(ymax,ynodes[edges[i,j]])
   end
   bbox[1]=xmin
   bbox[2]=xmax
   bbox[3]=ymin
   bbox[4]=ymax
end

function update_boundingbox!(bbox,xnodes,ynodes,edges,nnodes,edge_index)
   #updates bbox array with min(x),max(x),min(y),max(y)
   #only considers the cell in row edge_index of array edges 
   i=edge_index
   for j=1:nnodes[i]
      bbox[1]=min(bbox[1],xnodes[edges[i,j]])
      bbox[2]=max(bbox[2],xnodes[edges[i,j]])
      bbox[3]=min(bbox[3],ynodes[edges[i,j]])
      bbox[4]=max(bbox[4],ynodes[edges[i,j]])
   end
end

function update_boundingbox!(bbox1,bbox2)
   #updates bbox1 array with min(x),max(x),min(y),max(y)
   bbox1[1]=min(bbox1[1],bbox2[1])
   bbox1[2]=max(bbox1[2],bbox2[2])
   bbox1[3]=min(bbox1[3],bbox2[3])
   bbox1[4]=max(bbox1[4],bbox2[4])
end

function in_bbox(xpoint,ypoint,bbox::Array)
   return ( (xpoint>=bbox[1]) && (xpoint<=bbox[2]) && (ypoint>=bbox[3]) && (ypoint<=bbox[4]) )
end

#
# medium-level functions
#

function find_first_cell(xpoint,ypoint,grid::Grid)
   #for point (xpoint,ypoint) find cell with point insides
   #stops at first hit and return row number of edges array
   #return -1 if no cell is found
   #TODO start with linear search. Make more efficient later
   ##@printf("x=%f y=%f \n",xpoint,ypoint)
   if !in_bbox(xpoint,ypoint,grid.bbox)
      return -1 #no chance to find anything outside bbox
   end
   # first check subgrids for a match
   for igrid=1:length(grid.subgrids)
      icell=find_first_cell(xpoint,ypoint,grid.subgrids[igrid])
      if icell>0
         return icell
      end
   end
   #finally check local cells
   m=length(grid.edge_indices)
   wn=0
   s=0.
   x=0.
   for i=1:m
      iedge=grid.edge_indices[i]
      wn=0
      xprev=grid.xnodes[grid.edges[iedge,grid.nnodes[iedge]]]
      yprev=grid.ynodes[grid.edges[iedge,grid.nnodes[iedge]]]
      for j=1:grid.nnodes[iedge]
         ##@printf("edge %d\n",j)
         xthis=grid.xnodes[grid.edges[iedge,j]]
         ythis=grid.ynodes[grid.edges[iedge,j]]
         ##@printf("(%f,%f) -> (%f,%f)\n",xprev,yprev,xthis,ythis)
         # look for crossing of y=ypoint and x>xpoint, ie east of point
         if((yprev<ypoint)&&(ythis>=ypoint))
            #check for northward crossing east of point
   	    s=(ypoint-yprev)/(ythis-yprev)
   	    x=xprev+s*(xthis-xprev)
   	    if x>xpoint
   	       wn+=1
               ##@printf("north\n")
            end
         end
         if((yprev>=ypoint)&&(ythis<ypoint))
            #check for southward crossing east of point
   	    s=(ypoint-yprev)/(ythis-yprev)
   	    x=xprev+s*(xthis-xprev)
   	    if x>xpoint
   	       wn-=1
               ##@printf("south\n")
            end
         end
         xprev=xthis
         yprev=ythis
      end
      #wn now contains winding number
      ##@printf("cell %d winding-number %d \n",iedge,wn)
      if wn!=0 
         return iedge
      end
   end
   return -1
end


function find_cells(xpoints,ypoints,grid::Grid)
   #find cell numbers for each point in xpoints x ypoints
   #xpoints and ypoints are vectors with grid in 2 directions
   m=length(xpoints)
   n=length(ypoints)
   cells=zeros(Int64,m,n)
   for i=1:m
   for j=1:n
      cells[i,j]=find_first_cell(xpoints[i],ypoints[j],grid)
   end
   end
   return cells
end

function create_node_tree!(grid::Grid,nmin=100)
   if length(grid.edge_indices)>nmin
      # recursively split grid into subgrids
      n=length(grid.edge_indices)
      #@printf("n=%d\n",n)
      edge_indices1=zeros(Int64,n)
      edge_indices2=zeros(Int64,n)
      bbox1=zeros(Float64,4)
      bbox2=zeros(Float64,4)
      n1=0
      n2=0
      edge_indices=zeros(Int64,n)
      bbox=copy(grid.bbox)
      n=0
      #In which direction will we subdivide the box
      if( (bbox[2]-bbox[1]) > (bbox[4]-bbox[3]) )
         # split x-direction
         xsplit_low =bbox[1]+0.45*(bbox[2]-bbox[1])
         xsplit_high=bbox[1]+0.55*(bbox[2]-bbox[1])
         #@printf("x-split %f %f\n",xsplit_low,xsplit_high)
	 bbox_temp=zeros(Float64,4)
         for i=1:length(grid.edge_indices)
            iedge=grid.edge_indices[i]
            #@printf("edge %d\n",iedge)
            compute_boundingbox!(bbox_temp,grid.xnodes,grid.ynodes,grid.edges,grid.nnodes,iedge)
            if (bbox_temp[2]<=xsplit_high ) #left
               #println("left")
               if n1==0
                  bbox1=copy(bbox_temp)
               else
                  update_boundingbox!(bbox1,bbox_temp)
               end
               n1+=1
               edge_indices1[n1]=iedge
            elseif (bbox_temp[1] >= xsplit_low ) #right
               #println("right")
               if n2==0
                  bbox2=copy(bbox_temp)
               else
                  update_boundingbox!(bbox2,bbox_temp)
               end
               n2+=1
               edge_indices2[n2]=iedge
            else #keep at this level
               #println("between left and right")
               if n==0
                  bbox=copy(bbox_temp)
               else
                  update_boundingbox!(bbox,bbox_temp)
               end
               n+=1
               edge_indices[n]=iedge
            end
         end
      else
         # split y-direction
         ysplit_low =bbox[3]+0.45*(bbox[4]-bbox[3])
         ysplit_high=bbox[3]+0.55*(bbox[4]-bbox[3])
         #@printf("y-split %f %f\n",ysplit_low,ysplit_high)
	 bbox_temp=zeros(Float64,4)
         for i=1:length(grid.edge_indices)
            iedge=grid.edge_indices[i]
            #@printf("edge %d\n",iedge)
            compute_boundingbox!(bbox_temp,grid.xnodes,grid.ynodes,grid.edges,grid.nnodes,iedge)
            #println(bbox_temp)
            if (bbox_temp[4]<=ysplit_high ) #bottom
               #println("bottom")
               if n1==0
                  bbox1=copy(bbox_temp)
               else
                  update_boundingbox!(bbox1,bbox_temp)
               end
               n1+=1
               edge_indices1[n1]=iedge
            elseif (bbox_temp[3] >= ysplit_low ) #top
               #println("top")
               if n2==0
                  bbox2=copy(bbox_temp)
               else
                  update_boundingbox!(bbox2,bbox_temp)
               end
               n2+=1
               edge_indices2[n2]=iedge
            else #keep at this level
               #println("between top and bottom")
               if n==0
                  bbox=copy(bbox_temp)
               else
                  update_boundingbox!(bbox,bbox_temp)
               end
               n+=1
               edge_indices[n]=iedge
            end
         end
      end
      #@printf("parent,child1,child2=%d,%d,%d\n",n,n1,n2)
      #dump_array(bbox1,"bbox1")
      #dump_array(bbox2,"bbox2")
      #child 1
      resize!(edge_indices1,n1)
      #dump_array(edge_indices1,"edge_indices1")
      if(n1>0)
         grid1=Grid(grid.xnodes,grid.ynodes,grid.edges,edge_indices1,bbox1,grid.nnodes)
         create_node_tree!(grid1,nmin)
         push!(grid.subgrids,grid1)
      end
      #child 2
      resize!(edge_indices2,n2)
      #dump_array(edge_indices2,"edge_indices2")
      if(n2>0)
         grid2=Grid(grid.xnodes,grid.ynodes,grid.edges,edge_indices2,bbox2,grid.nnodes)
         create_node_tree!(grid2,nmin)
         push!(grid.subgrids,grid2)
      end
      #parent
      #dump_array(edge_indices,"edge_indices")
      resize!(edge_indices,n)
      grid.edge_indices=edge_indices
   end
end

import Base.dump

function dump(grid::Grid)
   dump_array(grid.xnodes,"xnodes")
   dump_array(grid.ynodes,"ynodes")
   @printf(" edges[] = ")
   dump(grid.edges)
   @printf("--- subgrid ---\n")
   @printf("%3d bbox=",0)
   show(grid.bbox)
   dump_array(grid.edge_indices,"  edge_indices")
   for i=1:length(grid.subgrids)
      dump_subgrid(grid.subgrids[i],1)
   end
end

function dump_subgrid(grid::Grid,level)
   @printf("%3d bbox=",level)
   show(grid.bbox)
   dump_array(grid.edge_indices,"  edge_indices")
   for i=1:length(grid.subgrids)
      dump_subgrid(grid.subgrids[i],level+1)
   end
end

function dump_array(a,name="")
   s=size(a)
   if(length(s)==1)
      @printf("%s[1:%d]=",name,length(a))
      if(length(a)>10)
         @printf("[%f,%f,%f,%f,%f,...,%f,%f,%f,%f,%f]\n",a[1],a[2],a[3],a[4],a[5],a[end-4],a[end-3],a[end-2],a[end-1],a[end])
      elseif(length(a)==0)
         print("[]\n")
      else
         show(a)
	 print("\n")
      end
   end
end

#
# testing 
# TODO move to proper test 
#

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
   return (xpoints0,ypoints0,wl_interp0)
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

#test1()
#test2()
#(x3,y3,z3)=test3()
#test4(20)

#(x5,y5,z5)=test5(10)
#using Plots
#contourf(5,y5,z5',clims=(-1.0,1.0))
#Profile.clear_malloc_data()
#@time test5(10)
#pcolormesh(xpoints,ypoints,wl_interp',vmin=-1,vmax=1)
