# Interact with dflow netcdf output map files
#

using NetCDF
using Dates
using Glob
# include("unstructured_grid.jl")

debuglevel = 1 # 0-nothing, larger more output


"""
   map=load_nc_info(path,filename)
Load meta-data of all map-files, i.e. for all domains in the filename.
"""
function load_nc_info(path, filename)
   map = []
   if isa(filename, Regex) # filename used to be a Regex-type (filename_regex::Regex)
      filenames = filter(x -> ~isnothing(match(filename, x)), readdir(path))
      filenames = [joinpath(path, filename) for filename = filenames]
   elseif occursin("*", filename)
      filenames = glob(filename, path)
   else
      filenames = [joinpath(path, filename)]
   end
   for filename = filenames
        @info filename
        push!(map, NetCDF.open(filename))
   end
   return map
end

"""
   interp=load_dflow_grid(map,nmin=50,spherical=true)
Create spatial index for a dflow netcdf map file.
map should contain an array of netcdf info's as created eg by:
load_nc_info.
"""
function load_dflow_grid(map, nmin=50, spherical=true)
   interp = Interpolator()
   println("compute index:")
   map_names = []
   for i = 1:length(map)
      # only load the grid with a certain filename once
      map_name = splitpath(map[i].name)[end]
      if in(map_name, map_names)
         break
      else
         push!(map_names, map_name)
      end
      if haskey(map[i].vars, "NetElemNode") # Old FM variable name
         edges_temp = map[i].vars["NetElemNode"][:,:]
      elseif haskey(map[i].vars, "mesh2d_face_nodes") # New FM variable name
	     edges_temp = map[i].vars["mesh2d_face_nodes"][:,:]
      elseif haskey(map[i].vars, "Mesh_face_nodes") # FEWS variable name
         edges_temp = map[i].vars["Mesh_face_nodes"][:,:]
      else
         error("Variable 'mesh2d_face_nodes' (or similar) is missing in D-Flow FM file")
      end
      if haskey(map[i].vars, "NetNode_x") # Old FM variable names
         xnodes_temp = map[i].vars["NetNode_x"][:]
         ynodes_temp = map[i].vars["NetNode_y"][:]
      elseif haskey(map[i].vars, "mesh2d_node_x") # New FM variable names
         xnodes_temp = map[i].vars["mesh2d_node_x"][:]
         ynodes_temp = map[i].vars["mesh2d_node_y"][:]
      elseif haskey(map[i].vars, "Mesh_node_x") # FEWS variable name
	      xnodes_temp = map[i].vars["Mesh_node_x"][:]
         ynodes_temp = map[i].vars["Mesh_node_y"][:]
      else
         error("Variable 'mesh2d_node_x' (or similar) is missing in D-Flow FM file")
      end
      println("- $(map[i].name)")
      grid_temp = Grid(xnodes_temp, ynodes_temp, edges_temp, nmin, spherical)
      # dump(grid_temp)
      add_grid!(interp, grid_temp)
   end
   return interp
end

"""
   depth = load_nc_var(map,"depth")
Load data of a variable foreach domain. The variable should be be time-independent.
"""
function load_nc_var(map, varname)
   result = []
   for i = 1:length(map)
      push!(result, map[i][varname][:])
   end
   return result
end

"""
   waterlevel = load_nc_map_slice(map,"s1",1)
Load data for a time-dependent variable for a specific time and for all domains.
"""
function load_nc_map_slice(map, varname, itime, sigma_layer_index=-1, domainnr = 0)
   result = []
   domainnr != 0 || (domainnr = 1:length(map))
   for i in domainnr
      # push!(result,map[i][varname][:,itime])
      if ndims(map[i][varname]) == 2
         push!(result, map[i][varname][:,itime])
      elseif ndims(map[i][varname]) == 3
         if sigma_layer_index > 0
            # if sigma_layer_index <= size(map[i][varname])[1]
               push!(result, map[i][varname][sigma_layer_index,:,itime])
            # else
            #    throw(BoundsError(map[i][varname], sigma_layer_index))
            # end
         else
            # For now only the upper sigma layer
            push!(result, map[i][varname][end,:,itime])
         end
      else
         throw(ArgumentError("load_nc_map_slice has only supports 2/3 dimensions"))
      end
   end
   return result
end


"""
   times=get_times(map,Dates.DateTime(2019,1,1))

Get available times im netcdf map files as a range in seconds, eg 0.0:3600.0:7200.0
The reftime (second arg) is a Dates.DateTime can be used to change the reftime to something
convenient. With reftime=get_reftime(map) you can find the time of the first map, which is often
a nice default.
"""
function get_times(map, reftime::DateTime)
   time_relative = map[1].vars["time"]
   units = time_relative.atts["units"]
   temp = split(units, "since")
   if endswith(temp[2],".0 +0000")
      t0 = DateTime(strip(temp[2]), "yyyy-mm-dd HH:MM:SS.0 +0000") # format used in FEWS 
   else
      t0 = DateTime(strip(temp[2]), "yyyy-mm-dd HH:MM:SS +00:00")
   end
   dt_seconds = 1.0
   if startswith(temp[1], "seconds")
      dt_seconds = 1.0
   elseif startswith(temp[1], "minutes")
      dt_seconds = 60.0
   elseif startswith(temp[1], "hours")
      dt_seconds = 3600.0
   elseif startswith(temp[1], "days")
      dt_seconds = 24.0 * 3600.0
   else
      error("Invalid time-step unit in map-file.")
   end
   # if ((time_relative[2]-time_relative[1])>(1.1*(time_relative[3]-time_relative[2])))
   #   #check for dummy initial field
   #   t_start = (0.001*Dates.value(t0-reftime))+time_relative[2]*dt_seconds
   #   t_step  = (time_relative[3]-time_relative[2])*dt_seconds
   #   t_stop  = (0.001*Dates.value(t0-reftime))+time_relative[end]*dt_seconds
   #   times = t_start:t_step:t_stop
   # else
   #   t_start = (0.001*Dates.value(t0-reftime))+time_relative[1]*dt_seconds
   #   t_step  = (time_relative[2]-time_relative[1])*dt_seconds
   #   t_stop  = (0.001*Dates.value(t0-reftime))+time_relative[end]*dt_seconds
   #   times = t_start:t_step:t_stop
   # end
   times = []
   for ti = 1:length(time_relative)
      push!(times, (0.001 * Dates.value(t0 - reftime)) + time_relative[ti] * dt_seconds)
   end
   return times
end

"""
   t_ref = get_reftime(map)
Read the reference time from the attributes of time in the netcdf file.
Times are described by a reference time (DateTime type) and the number of hours
relative to this t_ref.
Note that this function returns the first timestamp in the file (which 
   is not always the same as the RefDate from the .mdu or the 
   "since <refDate> text" in the time-attribute of the netcdf file)
Suggestion: use get_reftime to get the reference date and also use 
   get_reftime within get_times
"""
function get_reftime(map)
   time_relative = map[1].vars["time"]
   units = time_relative.atts["units"]
   temp = split(units, "since")
   if endswith(temp[2],"0 +00:00")
      t0 = DateTime(strip(temp[2]), "yyyy-mm-dd HH:MM:SS +00:00")
   else
      t0 = DateTime(strip(temp[2]), "yyyy-mm-dd HH:MM:SS")
   end
   dt_seconds = 1.0
   if startswith(temp[1], "seconds")
      dt_seconds = 1.0
   elseif startwith(temp[1], "minutes")
      dt_seconds = 60.0
   elseif startwith(temp[1], "hours")
      dt_seconds = 3600.0
   elseif startwith(temp[1], "days")
      dt_seconds = 24.0 * 3600.0
   else
      error("Invalid time-step unit in map-file.")
   end
   seconds_per_day = 24.0 * 3600.0
   # println("$(dt_seconds), $(t0), $(seconds_per_day), $(time_relative[1])")
   relative_days_in_seconds = div(dt_seconds * time_relative[1], seconds_per_day) * seconds_per_day
   return t0 + Dates.Second(relative_days_in_seconds)
end

"""
   t = as_DateTime(reftime,relative_time)
Convert a reftime (DateTime) and relative_time (Float64 with hours) relative to this.
Returns a DateTime type.
"""
function as_DateTime(reftime::DateTime, relative_time)
   reftime + Float64(relative_time) * Dates.Second(1)
end

"""
   u,v=initialize_interpolation(dflow_map,interp::Interpolator,varname,reftime::DateTime,dumval=0.0,cache_direction::Symbol=:forwards)

Create interpolation functions for x and y (lon and lat) directions of varname in dflow_map.
The reftime can be chosen freely independent of the reftime in the input files.
"""

function initialize_interpolation(dflow_map, interp::Interpolator, varname, reftime::DateTime, dumval=0.0, cache_direction::Symbol=:forwards)
   println("initialize caching for $(dflow_map[1].name) $varname...")
   # map_cache=dflow_map
   # interp_cache=interp
   # dumval_cache=dumval
   # reftime_cache=reftime
   times_cache = get_times(dflow_map, reftime)
   no_domains = length(interp.grids)
   multiple_runs = length(dflow_map) > no_domains # multiple runs per domain available
   if multiple_runs
      times_cache_per_run = [get_times([dflow_map[i]], reftime) for i in 1:no_domains:length(dflow_map)]
      dflow_map_all = dflow_map
      dflow_map_ind = 1:no_domains
      dflow_map = dflow_map[dflow_map_ind]
   end

   # keep 3 times in memory
   time_cache = zeros(3)
   time_cache_all_domains = zeros(3,no_domains)
   var_cache = Array{Any,1}(undef, 3)
   initialized = false
   if !haskey(dflow_map[1].vars, varname) 
      varname_stripped = replace(varname, "mesh2d_" => "") # e.g. mesh2d_ucx => ucx
      if haskey(dflow_map[1].vars, varname_stripped)
         println("Renaming varname to $varname_stripped as $varname is not a variable in the map, but $varname_stripped is.")
         varname = varname_stripped
      else
         throw(ArgumentError("$varname is not a variable in the map"))
      end
   end
   time_cache_all_domains = repeat(times_cache[1:3], 1, no_domains)
   time_cache_index_all_domains = repeat([3], no_domains) # index of last cached field in list of all available times
   for ti = 1:3
      var_cache[ti] = load_nc_map_slice(dflow_map, varname, ti)
   end
   (debuglevel > 4) && println("Initial cache index=$(time_cache_index_all_domains[1]) ")
   """
       update_cache_forwards(t)
       Refresh the cached fields in the forwards time directions if needed.
   """
   function update_cache_forwards(t, domainnr)
      time_cache = time_cache_all_domains[:, domainnr]
      time_cache_index = time_cache_index_all_domains[domainnr]

      # First, update dflow_map and refresh cache
      if multiple_runs
         run_ind_needed = findlast([(t >= times[1]) && (t <= times[end]) for times in times_cache_per_run])
         if !isnothing(run_ind_needed)
            dflow_map_ind_needed = ((run_ind_needed-1)*no_domains+1) : (run_ind_needed*no_domains)
            if dflow_map_ind_needed != dflow_map_ind # update of dflow_map needed (e.g. using run02/run_0000_map.nc instead of run01/run_0000_map.nc)
               (debuglevel >= 1) && println("t = $(t) - Update '$(varname)'-data using: ../$(joinpath(splitpath(dflow_map_all[dflow_map_ind_needed[domainnr]].name)[end-2:end]...))")
               dflow_map_ind = dflow_map_ind_needed
               dflow_map = dflow_map_all[dflow_map_ind]
               times_cache = times_cache_per_run[run_ind_needed]
               for ti = 1:3
                  time_cache[ti] = times_cache[ti]
                  var_cache[ti][domainnr] = load_nc_map_slice(dflow_map, varname, ti, -1, domainnr)[1]
               end
               time_cache_index = 3 # index of last cached field in list of all available times
            end
         end
      end

      # Next, update cache of this dflow_map 
      if (t >= time_cache[1]) && (t <= time_cache[3])
         (debuglevel >= 2) && println("cache is okay")
      elseif t > times_cache[end]
         error("Trying to access beyond last map t=$(t) > $(times_cache[end])")
      elseif (t >= time_cache[2]) && (t <= times_cache[time_cache_index + 1])
         (debuglevel >= 2) && println("advance to next time")
         time_cache[1] = time_cache[2]
         time_cache[2] = time_cache[3]
         time_cache[3] = times_cache[time_cache_index + 1]
         var_cache[1][domainnr] = var_cache[2][domainnr]
         var_cache[2][domainnr] = var_cache[3][domainnr]
         var_cache[3][domainnr] = load_nc_map_slice(dflow_map, varname, time_cache_index + 1, -1, domainnr)[1]
         time_cache_index += 1
      elseif t < times_cache[1]
         error("Trying to access before first map t=$(t) < $(times_cache[1])")
      else # complete refresh of cache
         (debuglevel >= 2) && println("refresh cache")
         ti = findfirst(tt -> tt > t, times_cache)
         (ti != length(times_cache)) || (ti = ti - 1)
         (ti != 1) || (ti = ti + 1)
         (debuglevel >= 4) && println("ti=$(ti), t=$(t)")
         (debuglevel >= 4) && println("$(times_cache)")
         time_cache[1] = times_cache[ti - 1]
         time_cache[2] = times_cache[ti]
         time_cache[3] = times_cache[ti + 1]
         var_cache[1][domainnr] = load_nc_map_slice(dflow_map, varname, ti - 1, -1, domainnr)[1]
         var_cache[2][domainnr] = load_nc_map_slice(dflow_map, varname, ti    , -1, domainnr)[1]
         var_cache[3][domainnr] = load_nc_map_slice(dflow_map, varname, ti + 1, -1, domainnr)[1]
         time_cache_index = ti + 1
      end
      (debuglevel >= 4) && println("$(time_cache_index) $(time_cache[1]) $(time_cache[2]) $(time_cache[3]) ")

      time_cache_all_domains[:, domainnr] = time_cache
      time_cache_index_all_domains[domainnr] = time_cache_index
   end

   """
       update_cache_backwards(t)
       Refresh the cached fields in the backwards time direction if needed.
   """
   function update_cache_backwards(t)
      if (t >= time_cache[1]) && (t <= time_cache[3])
         (debuglevel >= 2) && println("cache is okay")
      elseif t < times_cache[1]
         error("Trying to access before first map t=$(t) < $(times_cache[1])")
      elseif (t <= time_cache[2]) && (t >= times_cache[time_cache_index - 1])
         (debuglevel >= 2) && println("advance to next time")
         time_cache[3] = time_cache[2]
         time_cache[2] = time_cache[1]
         time_cache[1] = times_cache[time_cache_index - 1]
         var_cache[3] = var_cache[2]
         var_cache[2] = var_cache[1]
         var_cache[1] = load_nc_map_slice(dflow_map, varname, time_cache_index - 1)
         time_cache_index -= 1
      elseif t > times_cache[end]
         error("Trying to access beyond last map t=$(t) > $(times_cache[end])")
      else # complete refresh of cache
         (debuglevel >= 2) && println("refresh cache")
         ti = findfirst(tt -> tt > t, times_cache)
         (debuglevel >= 4) && println("ti=$(ti), t=$(t)")
         (debuglevel >= 4) && println("$(times_cache)")
         time_cache[1] = times_cache[ti - 2]
         time_cache[2] = times_cache[ti - 1]
         time_cache[3] = times_cache[ti]
         var_cache[1] = load_nc_map_slice(dflow_map, varname, ti - 2)
         var_cache[2] = load_nc_map_slice(dflow_map, varname, ti - 1)
         var_cache[3] = load_nc_map_slice(dflow_map, varname, ti)
         time_cache_index = ti - 2
      end
      (debuglevel >= 4) && println("$(time_cache_index) $(time_cache[1]) $(time_cache[2]) $(time_cache[3]) ")
      initialized = true
   end

   """
       (w1,w2,w3) = weights(t)
       Compute weights for (linear) time interpolation to time t based on 3 cached times
       This function assumes that the cache is up-to-date.
   """
   function weights(t)
      if t < time_cache[1] || t > time_cache[3]
         throw(ArgumentError("t outside cached time"))
      end
      if (t > time_cache[2])
         w = (t - time_cache[2]) / (time_cache[3] - time_cache[2])
         return (0.0, (1.0 - w), w)
      else
         w = (t - time_cache[1]) / (time_cache[2] - time_cache[1])
         return ((1.0 - w), w, 0.0)
      end
   end
   # flow in x direction (for now has to be called u)
   function f(x, y, z, t)
      ind = find_index(interp, x, y)
      if ind[1] == -1; return dumval; end
      if cache_direction == :forwards
         update_cache_forwards(t, ind[1])
      elseif cache_direction == :backwards
         update_cache_backwards(t)
      end
      w = weights(t)
      (debuglevel > 3) && println("weights $(weights)")
      value = 0.0
      for ti = 1:3
         temp = apply_index(ind, var_cache[ti], dumval)
         value += w[ti] * temp
         (debuglevel > 4) && println(" add w*val = $(w[ti]) $(temp)")
      end
      if abs(value) > 100.0 && debuglevel >= 2
         println("$varname interpolation")
         println("w $(weights) : $(ind)")
      end
      return value
   end
   return f
end

"""
   u,v=initialize_interpolation(dflow_map,interp::Interpolator,reftime::DateTime,dumval=0.0)

Create interpolation functions for x and y (lon and lat) directions of the flow.
The reftime can be chosen freely independent of the reftime in the input files.
"""
function initialize_interpolation(dflow_map, interp::Interpolator, reftime::DateTime, dumval=0.0)
   ucx = "ucx"
   ucy = "ucy"
   if haskey(dflow_map[1].vars, "mesh2d_ucx")
      ucx = "mesh2d_ucx"
      ucy = "mesh2d_ucy"
   end
   u = initialize_interpolation(dflow_map, interp, ucx, reftime, dumval)
   v = initialize_interpolation(dflow_map, interp, ucy, reftime, dumval)
   return u, v
end
