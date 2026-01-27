# some tools for unstructured grids
#
#CartesianGrid: spatial interpolation on regular x-y grid
# function CartesianGrid(xnodes::Array,ynodes::Array,spherical::Bool=true)
# function dump(grid::CartesianGrid)
# function in_bbox(grid::CartesianGrd,xpoint,ypoint)
# function find_index(grid::CartesianGrid,xpoint,ypoint)
# function find_index_and_weights(grid::CartesianGrid,xpoint,ypoint)
# function apply_index_and_weights(xindices,yindices,weights,values,dummy=0.0)
# function interpolate(grid::CartesianGrid,xpoint::Number,ypoint::Number,values,dummy=0.0)
# CartesianXYTGrid: space-time interpolation
# function CartesianXYTGrid(grid::CartesianGrid,times::AbstractVector,values::AbstractArrayi,name::String,missing_value::Number,scaling=1.0,offset=0.0)
# function interpolate(xyt::CartesianXYTGrid,xpoint::Number,ypoint::Number,time::Number,dummy=0.0)
# function get_map_slice(xyt::CartesianXYTGrid,ti:Integer)
# function update_cache(xyt::CartesianXYTGrid,t)
# function weights(xyt::CartesianXYTGrid,t)


debuglevel=1

if !@isdefined(CartesianGrid)
mutable struct CartesianGrid <: SpaceGrid
   xnodes::Array{Float64,1}     #x-coordinates of nodes
   ynodes::Array{Float64,1}     #y-coordinates of nodes
   regular::Bool # assume true for now
   #derived data
   bbox::Array{Float64,1}       #bounding box [xmin xmax ymin ymax]
   x0::Float64
   dx::Float64
   y0::Float64
   dy::Float64
   function CartesianGrid(xnodes::Array,ynodes::Array,spherical::Bool=true)
      bbox=zeros(Float64,4)
      bbox[1]=minimum(xnodes)
      bbox[2]=maximum(xnodes)
      bbox[3]=minimum(ynodes)
      bbox[4]=maximum(ynodes)
      x0=xnodes[1]
      dx=xnodes[2]-xnodes[1]
      y0=ynodes[1]
      dy=ynodes[2]-ynodes[1]
      return new(xnodes,ynodes,true,bbox,x0,dx,y0,dy)
   end
end
end #ifdef

#
# CartesianGrid functions
#
"""
   dump(grid)
Print a summary of the grid to screen.
"""
function dump(grid::CartesianGrid)
   println("bbox= $(grid.bbox)")
   print("xnodes=")
   display(grid.xnodes')
   print("ynodes=")
   display(grid.ynodes')
end

"""
   boolean_inbox=in_bbox(grid,xpoint,ypoint)
"""
function in_bbox(grid::CartesianGrid,xpoint,ypoint)
   result=( (xpoint>=grid.bbox[1]) && (xpoint<=grid.bbox[2])
	 && (ypoint>=grid.bbox[3]) && (ypoint<=grid.bbox[4]) )
   return result
end

"""
 m,n = find_index(grid,xpoint,ypoint)
Find index of point in the grid. The value is truncated downward.
"""
function find_index(grid::CartesianGrid,xpoint,ypoint)
   if in_bbox(grid,xpoint,ypoint)
      x_index=trunc(Int,((xpoint-grid.x0)/grid.dx))+1
      y_index=trunc(Int,((ypoint-grid.y0)/grid.dy))+1
      return (x_index,y_index)
   else
      return (-1,-1)
   end
end

"""
xindices,yindices,weights = find_index_and_weights(grid,xpoint,ypoint)
Get indices and weights for linear interpolation.
"""
function find_index_and_weights(grid::CartesianGrid,xpoint,ypoint)
   if in_bbox(grid,xpoint,ypoint)
      x_rel=(xpoint-grid.x0)/grid.dx
      y_rel=(ypoint-grid.y0)/grid.dy
      xi=trunc(Int,x_rel)+1
      yi=trunc(Int,y_rel)+1
      wx=x_rel-trunc(x_rel)
      wy=y_rel-trunc(y_rel)
      xindices=(xi,xi+1,xi+1,xi)
      yindices=(yi,yi,yi+1,yi+1)
      weights=((1-wx)*(1-wy),wx*(1-wy),wx*wy,(1-wx)*wy)
      return (xindices,yindices,weights)
   else
      return ((-1,-1,-1,-1),(-1,-1,-1,-1),(0.0,0.0,0.0,0.0))
   end
end

"""
value = apply_index_and_weights(xindiced,yindices,weights,values)
Compute interpolated value by applying the weights.
"""
function apply_index_and_weights(xindices,yindices,weights,values,dummy=0.0)
   if xindices[1]<0
      return dummy
   else
      result=0.0
      for i=1:length(xindices)
         result+=weights[i]*values[xindices[i],yindices[i]]
      end
      return result
   end
end

"""
"""
function interpolate(grid::CartesianGrid,xpoint::Number,ypoint::Number,values,dummy=0.0)
   xindices,yindices,weights = find_index_and_weights(grid,xpoint,ypoint)
   return apply_index_and_weights(xindices,yindices,weights,values,dummy)
end


#
# CartesianXYTGrid functions
#

#if !@isdefined(CartesianXYTGrid)
mutable struct CartesianXYTGrid <: SpaceTimeGrid
   grid::CartesianGrid  #grids for  domain
   times::AbstractVector
   values::AbstractArray  #source array for interpolation
   name::String
   #
   ndims::Int
   missing_value::Number #value used for non-valid values
   scaling::Float64
   offset::Float64
   # derived data
   cache::Array{Any,1}
   time_cache::Array{Float64,1}
   time_cache_index::Int64
   cache_direction::Symbol
   #constructor
   """
   xyt=CartesianXYTGrid(grid,times,values,"pressure",missing_value,scaling=1.0,offset=0.0)
   Create an xyt item for space-time interpolation.
   """
   function CartesianXYTGrid(grid::CartesianGrid,times::AbstractVector,values::AbstractArray,name::String,missing_value::Number,scaling=1.0,offset=0.0,cache_direction::Symbol=:forwards)
      (debuglevel>3) && println("initialize CartesianXYTGrid.")
      #keep 3 times in memmory
      time_cache=zeros(3)
      cache=Array{Any,1}(undef,3)
      ndims=length(size(values))-1 #time does not count
      for ti=1:3
         time_cache[ti]=times[ti]
         if ndims==2
            temp=values[:,:,ti]
         elseif ndims==3
            temp=values[:,:,:,ti]
         else
            error("Ndims should be 2 or 3 for a cartesian xyt-grids for now")
         end
         temp_scaled=offset.+scaling.*temp
         temp_scaled[temp.==missing_value].=NaN #use NaN for missing in cache
         cache[ti]=temp_scaled
      end
      time_cache_index=3 #index of last cached field in list of all available times
      (debuglevel>3) && println("Initial cache index=$(time_cache_index) ")
      if cache_direction!=:forwards&&cache_direction!=:backwards
         error("Unexpected symbol for cache_direction, $cache_direction not supported")
      end
      new(grid,times,values,name,ndims,missing_value,scaling,offset,cache,time_cache,time_cache_index,cache_direction)
   end
end
#end #ifdef

"""
slice = get_map_slice(xyt,ti)
Get time slice of dataset referred to in xyt for time index ti.
"""
function get_map_slice(xyt::CartesianXYTGrid,ti::Integer)
   if xyt.ndims==2
      temp=xyt.values[:,:,ti]
   elseif xyt.ndims==3
      temp=xyt.values[:,:,:,ti]
   end
   temp_scaled=xyt.offset.+xyt.scaling.*temp
   temp_scaled[temp.==xyt.missing_value].=NaN #use NaN for missing in cache
   return temp_scaled
end

"""
update_cache(xyt,t)
Advance cache to time t. Updating the content of xyt if necessary.
"""
function update_cache(xyt::CartesianXYTGrid,t)
   if (t>=xyt.time_cache[1])&&(t<=xyt.time_cache[3])
      (debuglevel>=2) && println("cache is okay")
   elseif (t>=xyt.time_cache[2])&&(t<=xyt.times[xyt.time_cache_index+1])
      if(t>xyt.times[end])
         error("Trying to access beyond last map t=$(t) > $(xyt.times_cache[end])")
      end
      (debuglevel>=2) && println("advance to next time")
      xyt.cache[1]=xyt.cache[2]
      xyt.cache[2]=xyt.cache[3]
      xyt.cache[3]=get_map_slice(xyt,xyt.time_cache_index+1)
      xyt.time_cache[1]=xyt.time_cache[2]
      xyt.time_cache[2]=xyt.time_cache[3]
      xyt.time_cache[3]=xyt.times[xyt.time_cache_index+1]
      xyt.time_cache_index+=1
   else #complete refresh of cache
      if(t>xyt.times[end])
         error("Trying to access beyond last map t=$(t) > $(xyt.times[end])")
      end
      if(t<xyt.times[1])
         error("Trying to access before first map t=$(t) < $(xyt.times[1])")
      end
      (debuglevel>=2) && println("refresh cache")
      ti=findfirst(tt->tt>t,xyt.times)
      (debuglevel>=4) && println("ti=$(ti), t=$(t)")
      (debuglevel>=4) && println("$(xyt.times)")
      xyt.time_cache[1]=xyt.times[ti-1]
      xyt.time_cache[2]=xyt.times[ti]
      xyt.time_cache[3]=xyt.times[ti+1]
      xyt.cache[1]=get_map_slice(xyt,ti-1)
      xyt.cache[2]=get_map_slice(xyt,ti)
      xyt.cache[3]=get_map_slice(xyt,ti+1)
      xyt.time_cache_index=ti+1
   end
   (debuglevel>=4) && println("$(xyt.time_cache_index) $(xyt.time_cache[1]) $(xyt.time_cache[2]) $(xyt.time_cache[3]) ")
end

"""
update_cache_backwards(xyt,t)
Advance cache to time t backwards in time. Updating the content of xyt if necessary.
"""
function update_cache_backwards(xyt::CartesianXYTGrid,t)
   if (t>=xyt.time_cache[1])&&(t<=xyt.time_cache[3])
      (debuglevel>=2) && println("cache is okay")
   elseif (t>=xyt.time_cache[2])&&(t<=xyt.times[xyt.time_cache_index+1])
      if(t>xyt.times[end])
         error("Trying to access beyond last map t=$(t) > $(xyt.times_cache[end])")
      end
      (debuglevel>=2) && println("advance to next time")
      xyt.cache[3]=xyt.cache[2]
      xyt.cache[2]=xyt.cache[1]
      xyt.cache[1]=get_map_slice(xyt,xyt.time_cache_index-1)
      xyt.time_cache[3]=xyt.time_cache[2]
      xyt.time_cache[2]=xyt.time_cache[1]
      xyt.time_cache[1]=xyt.times[xyt.time_cache_index-1]
      xyt.time_cache_index-=1
   else #complete refresh of cache
      if(t>xyt.times[end])
         error("Trying to access beyond last map t=$(t) > $(xyt.times[end])")
      end
      if(t<xyt.times[1])
         error("Trying to access before first map t=$(t) < $(xyt.times[1])")
      end
      (debuglevel>=2) && println("refresh cache")
      ti=findfirst(tt->tt>t,xyt.times)
      (debuglevel>=4) && println("ti=$(ti), t=$(t)")
      (debuglevel>=4) && println("$(xyt.times)")
      xyt.time_cache[1]=xyt.times[ti-2]
      xyt.time_cache[2]=xyt.times[ti-1]
      xyt.time_cache[3]=xyt.times[ti]
      xyt.cache[1]=get_map_slice(xyt,ti-2)
      xyt.cache[2]=get_map_slice(xyt,ti-1)
      xyt.cache[3]=get_map_slice(xyt,ti)
      xyt.time_cache_index=ti-2
   end
   (debuglevel>=4) && println("$(xyt.time_cache_index) $(xyt.time_cache[1]) $(xyt.time_cache[2]) $(xyt.time_cache[3]) ")
end

"""
(w1,w2,w3) = weights(xyt,t)
Compute weights for (linear) time interpolation to time t based on 3 cached times
This function assumes that the cache is up-to-date.
"""
function weights(xyt::CartesianXYTGrid,t)
	if t<xyt.time_cache[1]||t>xyt.time_cache[3]
		throw(ArgumentError("t outside cached time"))
	end
   if (t>xyt.time_cache[2])
      w=(t-xyt.time_cache[2])/(xyt.time_cache[3]-xyt.time_cache[2])
      return (0.0,(1.0-w),w)
   else
      w=(t-xyt.time_cache[1])/(xyt.time_cache[2]-xyt.time_cache[1])
      return ((1.0-w),w,0.0)
   end
end

"""
value = interpolate(xyt,x,y,t,0.0)
Interpolate in space and time.
"""
function interpolate(xyt::CartesianXYTGrid,xpoint::Number,ypoint::Number,time::Number,dummy=0.0)
   if xyt.cache_direction==:forwards
	    update_cache(xyt,time)
   elseif xyt.cache_direction==:backwards
    	update_cache_backwards(xyt,time)
   end
   w=weights(xyt,time)
   value=0.0
   if xyt.ndims==2
      for ti=1:3
         value+=w[ti]*interpolate(xyt.grid,xpoint,ypoint,xyt.cache[ti],NaN)
      end
   elseif xyt.ndims==3
      for ti=1:3
         value+=w[ti]*interpolate(xyt.grid,xpoint,ypoint,dropdims(xyt.cache[ti],dims=3),NaN)
      end
   end
   if isnan(value)
      value=dummy
   end
   return value
end
