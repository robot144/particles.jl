# some tools for unstructured grids
using Test
using NetCDF

debuglevel=0

if !@isdefined(Grid)
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
end #ifdef

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

"""
function interpolate(interp::Interpolator,xpoints,ypoints,invalues::Array,dumval=0.0)
   interpolation for a set of points. The indices of the points are cached if they stay the say between calls
"""
function interpolate(interp::Interpolator,xpoints,ypoints,invalues::Array,dumval=0)
   #check for mutiple or single domain
   if(typeof(invalues).parameters[1] <: Number)
      #input as single domain
      if length(interp.xnodes)>1
         error("Input values were given for single domain to multi-domain interpolator")
      end
      values=[invalues]
   else
      values=invalues
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
         println("grid $(i)")
         @time cells=find_cells(xpoints,ypoints,interp.grids[i])
         push!(cached_indices,cells)
      end
      interp.cached_indices=cached_indices #store to cache for next call to this function
   else
      if debuglevel>1
         println("use cached interpolation weights")
      end
      cached_indices=interp.cached_indices #retrieve from cache
   end
   outvalues=zeros(Float64,size(cached_indices[1]))
   outvalues.=dumval
   for i=1:length(cached_indices) #loop over grids
      get_values_by_cells!(outvalues,cached_indices[i],values[i])
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
   ##println("x=$(xpoint) y=$(ypoint) ")
   m,n=size(edges)
   for i=1:m
      wn=0
      xprev=xnodes[edges[i,nnodes[i]]]
      yprev=ynodes[edges[i,nnodes[i]]]
      for j=1:nnodes[i]
         ##println("edge $(j)")
         xthis=xnodes[edges[i,j]]
         ythis=ynodes[edges[i,j]]
         ##println("($(xprev),$(yprev))) -> ($(xthis),$(ythis))")
         # look for crossing of y=ypoint and x>xpoint, ie east of point
         if((yprev<ypoint)&&(ythis>=ypoint))
            #check for northward crossing east of point
            s=(ypoint-yprev)/(ythis-yprev)
   	    x=xprev+s*(xthis-xprev)
   	    if x>xpoint
   	       wn+=1
               ##println("north")
            end
         end
         if((yprev>=ypoint)&&(ythis<ypoint))
            #check for southward crossing east of point
   	    s=(ypoint-yprev)/(ythis-yprev)
   	    x=xprev+s*(xthis-xprev)
   	    if x>xpoint
   	       wn-=1
               ##println("south")
            end
         end
         xprev=xthis
         yprev=ythis
      end
      #wn now contains winding number
      ##println("cell $(i) winding-number $(wn) ")
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
   ##println("x=$(xpoint) y=$(ypoint)")
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
         ##println("edge $(i))")
         xthis=grid.xnodes[grid.edges[iedge,j]]
         ythis=grid.ynodes[grid.edges[iedge,j]]
         ##println("($(xprev),$(yprev)) -> ($(xthis),$(ythis))\n",xprev,yprev,xthis,ythis)
         # look for crossing of y=ypoint and x>xpoint, ie east of point
         if((yprev<ypoint)&&(ythis>=ypoint))
            #check for northward crossing east of point
   	    s=(ypoint-yprev)/(ythis-yprev)
   	    x=xprev+s*(xthis-xprev)
   	    if x>xpoint
   	       wn+=1
               ##println("north")
            end
         end
         if((yprev>=ypoint)&&(ythis<ypoint))
            #check for southward crossing east of point
   	    s=(ypoint-yprev)/(ythis-yprev)
   	    x=xprev+s*(xthis-xprev)
   	    if x>xpoint
   	       wn-=1
               ##println("south")
            end
         end
         xprev=xthis
         yprev=ythis
      end
      #wn now contains winding number
      ##println("cell $(iednge) winding-number $(wn)")
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
      #println("n=$(n)")
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
         #println("x-split $(xsplit_low) $(xsplit_high)")
	 bbox_temp=zeros(Float64,4)
         for i=1:length(grid.edge_indices)
            iedge=grid.edge_indices[i]
            #println("edge $(iedge)")
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
         #println("y-split $(ysplit_low) $(ysplit_high)")
	 bbox_temp=zeros(Float64,4)
         for i=1:length(grid.edge_indices)
            iedge=grid.edge_indices[i]
            #print("edge $(iedge)")
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
      #println("parent,child1,child2=$(n),$(n1),$(n2)")
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
   print(" edges[] = ")
   dump(grid.edges)
   println("--- subgrid ---")
   print("000 bbox=")
   show(grid.bbox)
   dump_array(grid.edge_indices,"  edge_indices")
   for i=1:length(grid.subgrids)
      dump_subgrid(grid.subgrids[i],1)
   end
end

function dump_subgrid(grid::Grid,level)
   print("$(level) bbox=")
   show(grid.bbox)
   dump_array(grid.edge_indices,"  edge_indices")
   for i=1:length(grid.subgrids)
      dump_subgrid(grid.subgrids[i],level+1)
   end
end

function dump_array(a,name="")
   s=size(a)
   if(length(s)==1)
      print("$(name)[1:$(length(a))]=")
      if(length(a)>10)
         println("[$(a[1]),$(a[2]),$(a[3]),$(a[4]),$(a[5]),...,
		   $(a[end-4]),$(a[end-3]),$(a[end-2]),$(a[end-1]),$(a[end])]")
      elseif(length(a)==0)
         println("[]")
      else
         show(a)
	 println("")
      end
   end
end

"""
   cell_index = find_index(interp, xpoint, ypoint)
Find the domain and index of the cell within that domain, eg the result [2,1234]
indicates the cell 1234 in the snd domain.
"""
function find_index(interp::Interpolator,xpoint,ypoint)
   indices=[-1 -1]
   for i=1:length(interp.grids)
      cell=find_first_cell(xpoint,ypoint,interp.grids[i])
      if cell>0
         indices[1]=i
         indices[2]=cell
         break
      end
   end
   return indices
end

"""
   waterlevel_at_point = apply_index(index,map_slice,9999.0)
Get the number at domain and index (given by index). The index is often the result of 
the function find_index. If the cell index is [-1,-1] then the point is considered to
be outside the area covered by the cells, eg on land, and then a default value is returned.
"""
function apply_index(index,map_slice,default_value=0.0)
   if index[1]>0
      return map_slice[index[1]][index[2]]
   else
      return default_value
   end
end
