# Interact with dflow netcdf output map files
#

using NetCDF
using Dates
#include("unstructured_grid.jl")

const debuglevel=1 #0-nothing, larger more output


"""
   map=load_nc_info(path,filename_regex)
Load meta-data of all map-files, i.e. for all domains in the filename_regex.
"""
function load_nc_info(path,filename_regex::Regex)
   map=[]
   filenames=filter(x->(match(filename_regex,x)!=nothing), readdir(path))
   for filename=filenames
      push!(map,NetCDF.open(joinpath(path,filename)))
   end
   return map
end

"""
   interp=load_dflow_grid(map,nmin=50,spherical=true)
Create spatial index for a dflow netcdf map file.
map should contain an array of netcdf info's as created eg by:
load_nc_info.
"""
function load_dflow_grid(map,nmin=50,spherical=true)
   interp=Interpolator()
   println("compute index:")
   for i=1:length(map)
      if haskey(map[i].vars,"NetElemNode")
         edges_temp=map[i].vars["NetElemNode"][:,:]
      else
         edges_temp=map[i].vars["mesh2d_face_nodes"][:,:]
      end
      if haskey(map[i].vars,"NetNode_x")
         xnodes_temp=map[i].vars["NetNode_x"][:]
         ynodes_temp=map[i].vars["NetNode_y"][:]
      else
         xnodes_temp=map[i].vars["mesh2d_node_x"][:]
         ynodes_temp=map[i].vars["mesh2d_node_y"][:]
      end
      println("- $(map[i].name)")
      grid_temp=Grid(xnodes_temp,ynodes_temp,edges_temp,nmin,spherical)
      #dump(grid_temp)
      add_grid!(interp,grid_temp)
   end
   return interp
end

"""
   depth = load_nc_var(map,"depth")
Load data of a variable foreach domain. The variable should be be time-independent.
"""
function load_nc_var(map,varname)
   result=[]
   for i=1:length(map)
      push!(result,map[i][varname][:])
   end
   return result
end

"""
   waterlevel = load_nc_map_slice(map,"s1",1)
Load data for a time-dependent variable for a specific time and for all domains.
"""
function load_nc_map_slice(map,varname,itime)
   result=[]
   for i=1:length(map)
      push!(result,map[i][varname][:,itime])
   end
   return result
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

"""
   times=get_times(map,Dates.DateTime(2019,1,1))

Get available times im netcdf map files as a range in seconds, eg 0.0:3600.0:7200.0
The reftime (second arg) is a Dates.DateTime can be used to change the reftime to something
convenient. With reftime=get_reftime(map) you can find the time of the first map, which is often
a nice default.
"""
function get_times(map,reftime::DateTime)
   time_relative=map[1].vars["time"]
   units=time_relative.atts["units"]
   temp=split(units,"since")
   t0=DateTime(strip(temp[2]),"yyyy-mm-dd HH:MM:SS")
   dt_seconds=1.0
   if startswith(temp[1],"seconds")
      dt_seconds=1.0
   elseif startwith(temp[1],"minutes")
      dt_seconds=60.0
   elseif startwith(temp[1],"hours")
      dt_seconds=3600.0
   elseif startwith(temp[1],"days")
      dt_seconds=24.0*3600.0
   else
      error("Invalid time-step unit in map-file.")
   end
   if ((time_relative[2]-time_relative[1])>(1.1*(time_relative[3]-time_relative[2])))
      #check for dummy initial field
      t_start = (0.001*Dates.value(t0-reftime))+time_relative[2]*dt_seconds
      t_step  = (time_relative[3]-time_relative[2])*dt_seconds
      t_stop  = (0.001*Dates.value(t0-reftime))+time_relative[end]*dt_seconds
      times = t_start:t_step:t_stop
   else
      t_start = (0.001*Dates.value(t0-reftime))+time_relative[1]*dt_seconds
      t_step  = (time_relative[2]-time_relative[1])*dt_seconds
      t_stop  = (0.001*Dates.value(t0-reftime))+time_relative[end]*dt_seconds
      times = t_start:t_step:t_stop
   end
   return times
end

"""
   t_ref = get_reftime(map)
Read the reference time from the attributes of time in the netcdf file. 
Times are described by a reference time (DateTime type) and the number of hours
relative to this t_ref.
"""
function get_reftime(map)
   time_relative=map[1].vars["time"]
   units=time_relative.atts["units"]
   temp=split(units,"since")
   t0=DateTime(strip(temp[2]),"yyyy-mm-dd HH:MM:SS")
   dt_seconds=1.0
   if startswith(temp[1],"seconds")
      dt_seconds=1.0
   elseif startwith(temp[1],"minutes")
      dt_seconds=60.0
   elseif startwith(temp[1],"hours")
      dt_seconds=3600.0
   elseif startwith(temp[1],"days")
      dt_seconds=24.0*3600.0
   else
      error("Invalid time-step unit in map-file.")
   end
   seconds_per_day=24.0*3600.0
   #println("$(dt_seconds), $(t0), $(seconds_per_day), $(time_relative[1])")
   relative_days_in_seconds=div(dt_seconds*time_relative[1],seconds_per_day)*seconds_per_day
   return t0 + Dates.Second(relative_days_in_seconds)
end

"""
   t = as_DateTime(reftime,relative_time)
Convert a reftime (DateTime) and relative_time (Float64 with hours) relative to this.
Returns a DateTime type.
"""
function as_DateTime(reftime::DateTime,relative_time)
   reftime+Float64(relative_time)*Dates.Second(1)
end

"""
   u,v=initialize_interpolation(dflow_map,interp::Interpolator,reftime::DateTime,dumval=-9999.0)

Create interpolation functions for x and y (lon and lat) directions of the flow. 
The reftime can be chosen freely independent of the reftime in the input files.
""" 
function initialize_interpolation(dflow_map,interp::Interpolator,reftime::DateTime,dumval=0.0)
   println("initialize caching for $(dflow_map[1].name) ...")
   map_cache=dflow_map
   interp_cache=interp
   dumval_cache=dumval
   reftime_cache=reftime
   times_cache=get_times(map_cache,reftime_cache)
   #keep 3 times in memmory
   time_cache=zeros(3)
   u_cache=Array{Any,1}(undef,3)
   v_cache=Array{Any,1}(undef,3)
   initialized=false
   ucx="ucx"
   ucy="ucy"
   if haskey(dflow_map[1].vars,"mesh2d_ucx")
      ucx="mesh2d_ucx"
      ucy="mesh2d_ucy"
   end
   for ti=1:3
      time_cache[ti]=times_cache[ti]
      u_cache[ti]=load_nc_map_slice(dflow_map,ucx,ti)
      v_cache[ti]=load_nc_map_slice(dflow_map,ucy,ti)
   end
   time_cache_index=3 #index of last cached field in list of all available times
   """
       update_cache(t)
       Refresh the cached fields if needed.
   """
   function update_cache(t)
      ucx="ucx"
      ucy="ucy"
      if haskey(dflow_map[1].vars,"mesh2d_ucx")
         ucx="mesh2d_ucx"
         ucy="mesh2d_ucy"
      end
      if (t>=time_cache[1])&&(t<=time_cache[3])
         (debuglevel>=2) && println("cache is okay")
      elseif (t>=time_cache[2])&&(t<=times_cache[time_cache_index+1]) 
         if(t>times_cache[end])
            error("Trying to access beyond last map t=$(t) > $(times_cache[end])")
         end
         (debuglevel>=2) && println("advance to next time")
         u_cache[1]=u_cache[2]
         u_cache[2]=u_cache[3]
         u_cache[3]=load_nc_map_slice(dflow_map,ucx,time_cache_index+1)
         v_cache[1]=v_cache[2]
         v_cache[2]=v_cache[3]
         v_cache[3]=load_nc_map_slice(dflow_map,ucy,time_cache_index+1)
         time_cache[1]=time_cache[2]
         time_cache[2]=time_cache[3]
         time_cache[3]=times_cache[time_cache_index+1]
         time_cache_index+=1
      else #complete refresh of cache
         if(t>times_cache[end])
            error("Trying to access beyond last map t=$(t) > $(times_cache[end])")
         end
         if(t<times_cache[1])
            error("Trying to access before first map t=$(t) < $(times_cache[0])")
         end
         (debuglevel>=2) && println("refresh cache")
         ti=findfirst(tt->tt>=t,times_cache)
         (debuglevel>=4) && println("ti=$(ti), t=$(t)")
         (debuglevel>=4) && println("$(times_cache)")
         time_cache[1]=times_cache[ti]
         time_cache[2]=times_cache[ti+1]
         time_cache[3]=times_cache[ti+2]
         u_cache[1]=load_nc_map_slice(dflow_map,ucx,ti)
         u_cache[2]=load_nc_map_slice(dflow_map,ucx,ti+1)
         u_cache[3]=load_nc_map_slice(dflow_map,ucx,ti+2)
         v_cache[1]=load_nc_map_slice(dflow_map,ucy,ti)
         v_cache[2]=load_nc_map_slice(dflow_map,ucy,ti+1)
         v_cache[3]=load_nc_map_slice(dflow_map,ucy,ti+2)
         time_cache_index=ti+2
      end
      (debuglevel>=4) && println("$(time_cache_index) $(time_cache[1]) $(time_cache[2]) $(time_cache[3]) ")
   end
   """
       (w1,w2,w3) = weights(t)
       Compute weights for (linear) time interpolation to time t based on 3 cached times
       This function assumes that the cache is up-to-date.
   """
   function weights(t)
      if (t>time_cache[2])
         w=(t-time_cache[2])/(time_cache[3]-time_cache[2])
         return (0.0,(1.0-w),w)
      else
         w=(t-time_cache[1])/(time_cache[2]-time_cache[1])
         return ((1.0-w),w,0.0)
      end
   end
   #flow in x direction (for now has to be called u)
   function u(x,y,z,t)
      ind=find_index(interp_cache,x,y)
      update_cache(t)
      w=weights(t)
      value=0.0
      for ti=1:3
         value+=w[ti]*apply_index(ind,u_cache[ti],dumval_cache)
      end
      if abs(value)>100.0
         println("u interpolation")
         println("w $(weights) : $(ind)")
      end
      return value
   end
   #flow in y direction (for now has to be called v)
   function v(x,y,z,t)
      ind=find_index(interp_cache,x,y)
      update_cache(t)
      w=weights(t)
      value=0.0
      for ti=1:3
         value+=w[ti]*apply_index(ind,v_cache[ti],dumval_cache)
      end
      return value
   end
   return (u,v)
end


