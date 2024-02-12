# Interact with ERA5 netcdf output
# basically these are lon-lat regular NetCDF files.
# These routines are tested with data from the copernicus
# climate data store
# For more see: https://cds.climate.copernicus.eu
#
# function EraData(path,filename)
# function load_map_slice(data::EraData,varname,itime)
# function get_reftime(data::EraData)
# function get_times(data::EraData,reftime::DateTime)
# function as_DateTime(data::EraData,reftime::DateTime,relative_time)
# function initialize_interpolation(data::EraData,varname::String,reftime::DateTime,dummy=0.0)


using NetCDF
using Dates

debuglevel=1 #0-nothing, larger more output

"""
Using the EraData struct the file-handling becomes object-oriented.
"""
mutable struct EraData
   file::NcFile
   #derived data
   grid::CartesianGrid
   """
   Constructor
   era_data = EraData(".","my_era_file.nc")
   """
   function EraData(path,filename)
      file=NetCDF.open(joinpath(path,filename))
      x=collect(file.vars["longitude"])
      y=collect(file.vars["latitude"])
      grid=CartesianGrid(x,y,true)
      return new(file,grid)
   end
end

"""
   wind_x = load_map_slice(era_data,"u10",1)
Load data for a time-dependent variable for a specific time.
"""
function load_map_slice(data::EraData,varname,itime)
   ndims=length(size(data.file[varname]))-1 #only spatial dims
   scale_factor=data.file[varname].atts["scale_factor"]
   offset=data.file[varname].atts["add_offset"]
   dummy=data.file[varname].atts["missing_value"]
   if ndims==2
      tempvar=data.file[varname][:,:,itime]
      if length(size(tempvar))>2 # Starting somewhere around julia version 1.3 dimensions are dropped automatically when slicing
         # this check can be removed some time in the future, when we assumen versions less than 1.4 are no longer in use
         tempvar=dropdims(tempvar,dims=3)
      end
      var=offset.+scale_factor.*tempvar
      var[tempvar.==dummy].=NaN
      return var
   else
      error("only 2D (one layer) variables supported for now.")
   end
end

"""
   t_ref = get_reftime(era_data)
Read the reference time from the attributes of time in the netcdf file.
Times are described by a reference time (DateTime type) and the number of hours
relative to this t_ref.
"""
function get_reftime(data::EraData)
   time_relative=data.file.vars["time"]
   units=time_relative.atts["units"]
   temp=split(units,"since")
   temp2=split(temp[2],".")
   println("ref date $(temp2[1])")
   t0=DateTime(strip(temp2[1]),"yyyy-mm-dd HH:MM:SS")
   dt_seconds=1.0
   if startswith(temp[1],"seconds")
      dt_seconds=1.0
   elseif startswith(temp[1],"minutes")
      dt_seconds=60.0
   elseif startswith(temp[1],"hours")
      dt_seconds=3600.0
   elseif startswith(temp[1],"days")
      dt_seconds=24.0*3600.0
   else
      error("Invalid time-step unit in map-file.")
   end
   seconds_per_day=24.0*3600.0
   if debuglevel>2
      println("$(dt_seconds), $(t0), $(seconds_per_day), $(time_relative[1])")
   end
   relative_days_in_seconds=div(dt_seconds*time_relative[1],seconds_per_day)*seconds_per_day
   return t0 + Dates.Second(relative_days_in_seconds)
end

"""
   times=get_times(era_data,Dates.DateTime(2019,1,1))

Get available times im netcdf map files as a range in seconds, eg 0.0:3600.0:7200.0
The reftime (second arg) is a Dates.DateTime can be used to change the reftime to something
convenient. With reftime=get_reftime(map) you can find the time of the first map, which is often
a nice default.
"""
function get_times(data::EraData,reftime::DateTime)
   time_relative=data.file.vars["time"]
   units=time_relative.atts["units"]
   temp=split(units,"since")
   temp2=split(temp[2],".")
   println("ref date $(temp2[1])")
   t0=DateTime(strip(temp2[1]),"yyyy-mm-dd HH:MM:SS")
   dt_seconds=1.0
   if startswith(temp[1],"seconds")
      dt_seconds=1.0
   elseif startswith(temp[1],"minutes")
      dt_seconds=60.0
   elseif startswith(temp[1],"hours")
      dt_seconds=3600.0
   elseif startswith(temp[1],"days")
      dt_seconds=24.0*3600.0
   else
      error("Invalid time-step unit in map-file.")
   end
   times=[]
   for ti=1:length(time_relative)
      push!(times,(0.001*Dates.value(t0-reftime))+time_relative[ti]*dt_seconds)
   end
   return times
end

"""
   t = as_DateTime(data,reftime,relative_time)
Convert a reftime (DateTime) and relative_time (Float64 with seconds) relative to this.
Returns a DateTime type.
"""
function as_DateTime(data::EraData,reftime::DateTime,relative_time)
   reftime+Float64(relative_time)*Dates.Second(1)
end

"""
p = initialize_interpolation(era,"msl",t0)
Create an interpolation function p(x,y,z,t)
"""
function initialize_interpolation(data::EraData,varname::String,reftime::DateTime,dummy=0.0,cache_direction::Symbol=:forwards)
   times=get_times(data,reftime)
   values=data.file.vars[varname] #TODO more checks
   missing_value=values.atts["missing_value"]
   scaling=values.atts["scale_factor"]
   offset=values.atts["add_offset"]
   xyt=CartesianXYTGrid(data.grid,times,values,varname,missing_value,scaling,offset,cache_direction)
   function f(x,y,z,t)
      value=interpolate(xyt,x,y,t,dummy)
      return value
   end
   return f
end
