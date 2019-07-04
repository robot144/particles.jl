# Interact with dflow files
#


using NetCDF
using Dates
include("unstructured_grid.jl")


function load_nc_info(path,filename_regex::Regex)
   map=[]
   filenames=filter(x->(match(filename_regex,x)!=nothing), readdir(path))
   for filename=filenames
      push!(map,ncinfo(joinpath(path,filename)))
   end
   return map
end

function load_dflow_grid(map,nmin=50,spherical=true)
   interp=Interpolator()
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
      @printf("- index computation\n")
      @time grid_temp=Grid(xnodes_temp,ynodes_temp,edges_temp,nmin,spherical)
      #dump(grid_temp)
      add_grid!(interp,grid_temp)
   end
   return interp
end

function load_nc_var(map,varname)
   result=[]
   for i=1:length(map)
      push!(result,map[i][varname][:])
   end
   return result
end

function load_nc_map_slice(map,varname,itime)
   result=[]
   for i=1:length(map)
      push!(result,map[i][varname][:,itime])
   end
   return result
end

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

function apply_index(index,map_slice,default_value=0.0)
   if index[1]>0
      return map_slice[index[1]][index[2]]
   else
      return default_value
   end
end

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
      t_start = (0.001*Dates.value(reftime-t0))+time_relative[2]*dt_seconds
      t_step  = (time_relative[3]-time_relative[2])*dt_seconds
      t_stop  = (0.001*Dates.value(reftime-t0))+time_relative[end]*dt_seconds
      times = t_start:t_step:t_stop
   else
      t_start = (0.001*Dates.value(reftime-t0))+time_relative[1]*dt_seconds
      t_step  = (time_relative[2]-time_relative[1])*dt_seconds
      t_stop  = (0.001*Dates.value(reftime-t0))+time_relative[end]*dt_seconds
      times = t_start:t_step:t_stop
   end
   return times
end

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

function as_DateTime(reftime::DateTime,relative_time)
   reftime+Float64(relative_time)*Dates.Second(1)
end

function initialize_interpolation(dflow_map,interp::Interpolator,reftime::DateTime,dumval=-9999.0)
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
         println("cache is okay")
      elseif (t>=time_cache[2])&&(t<=times_cache[time_cache_index+1]) 
         println("advance to next time")
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
         println("refresh cache")
         ti=findfirst(tt->tt>=t,times_cache)
         println("ti=$(ti)")
         println("$(times_cache)")
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
      println("$(time_cache_index) $(time_cache[1]) $(time_cache[2]) $(time_cache[3]) ")
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
         value+=apply_index(ind,u_cache[ti],dumval_cache)
      end
      value
   end
   #flow in y direction (for now has to be called v)
   function v(x,y,z,t)
      ind=find_index(interp_cache,x,y)
      update_cache(t)
      w=weights(t)
      value=0.0
      for ti=1:3
         value+=apply_index(ind,v_cache[ti],dumval_cache)
      end
      value
   end
   return (u,v)
end
#
# test
#

function test1()
   #init
   dflow_map=load_nc_info("../test_data",r"estuary_...._map.nc")
   interp=load_dflow_grid(dflow_map,50,false)
   
   #interpolate a field to a regular grid
   sealevel=load_nc_map_slice(dflow_map,"s1",10)
   u_velocity=load_nc_map_slice(dflow_map,"ucx",10)
   xpoints=collect(range(0.,stop=100000.,length=300))
   ypoints=collect(range(0.,stop=500.,length=100))
   sealevel_interp=interpolate(interp,xpoints,ypoints,sealevel)
   u_velocity_interp=interpolate(interp,xpoints,ypoints,u_velocity)
   #heatmap(xpoints,ypoints,u_velocity_interp')
   # interpolate for one point only
   ind=find_index(interp,10000.0,200.0)
   ux=apply_index(ind,u_velocity,-9999.)
   
   # u,v interpolation functions
   t0=get_reftime(dflow_map)
   u,v=initialize_interpolation(dflow_map,interp,t0)
   u(100.0,100.0,0.0,0.0)
end

function test2()
   #init
   dflow_map=load_nc_info("../test_data",r"DCSM-FM_0_5nm_...._map.nc")
   interp=load_dflow_grid(dflow_map,50,false)
   ##interpolate a field to a regular grid
   #sealevel=load_nc_map_slice(dflow_map,"mesh2d_s1",10)
   #xpoints=collect(range(-15.0,stop=13.0,length=1200))
   #ypoints=collect(range(43.0,stop=64.0,length=1000))
   #sealevel_interp=interpolate(interp,xpoints,ypoints,sealevel)
   #Plots.default(:size,[1200,1000])
   #heatmap(xpoints,ypoints,sealevel_interp', clims=(-2.0,2.0))

   # u,v interpolation functions
   t0=get_reftime(dflow_map)
   u,v=initialize_interpolation(dflow_map,interp,t0)
   for istep=1:5
      u_value=u(1.0,51.0,0.0,864000.0+3600.0*istep)
      println("$(istep) u=$(u_value)")
   end
end

#test1()
#test2()
