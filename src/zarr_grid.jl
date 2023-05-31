# IntZarrct with gridded data in Zarr files
# basically these are regulat lon-lat or x-y files.
# The zarr format allows for chunking and compression and can be accessed
# on the local filesystem or via http.
# Note that the Zarr format uses multiple files in a directory structure, with one chunk per file.
#
# function ZarrData(path,filename)
# function load_map_slice(data::ZarrData,varname,itime)
# function get_reftime(data::ZarrData)
# function get_times(data::ZarrData,reftime::DateTime)
# function as_DateTime(data::ZarrData,reftime::DateTime,relative_time)
# function initialize_interpolation(data::ZarrData,varname::String,reftime::DateTime,dummy=0.0)


using Zarr
using Dates
using CFTime

debuglevel=1 #0-nothing, larger more output

"""
Using the ZarrData struct the file-handling becomes object-oriented.
"""
mutable struct ZarrData
   file::ZGroup
   #derived data
   grid::CartesianGrid
   """
   Constructor
   Zarr_data = ZarrData(".","my_Zarr_file.zarr")
   """
   function ZarrData(path,filename,spherical=true)
      file=zopen(joinpath(path,filename))
      x=collect(file.arrays["x_center"])
      y=collect(file.arrays["y_center"])
      grid=CartesianGrid(x,y,spherical)
      return new(file,grid)
   end
end

"""
    Zarr_Data = ZarrData(".",filename)
    names = varnames(Zarr_data)
    "time" in names
"""
function varnames(data::ZarrData)
   return keys(data.file.arrays)
end

"""
    z=ZarrData(".","my_Zarr_file.zarr")
    "time" in get_varnames(z)
"""
function get_varnames(data::ZarrData)
   return keys(data.file.arrays)
end


"""
   u = load_map_slice(Zarr_data,"x_velocity",1)
Load data for a time-dependent variable for a specific time.
"""
function load_map_slice(data::ZarrData,varname,itime)
    if !(varname in get_varnames(data))
        error("variable $varname not found in $(get_varnames(data))")
    end
    ndims=length(size(data.file[varname]))-1 #only spatial dims, no time
    scale_factor=data.file[varname].attrs["scale_factor"]
    offset=data.file[varname].attrs["add_offset"]
    dummy=data.file[varname].attrs["_FillValue"]
    if ndims==2
        tempvar=data.file[varname][:,:,itime]
        var=offset.+scale_factor.*tempvar
        var[tempvar.==dummy].=NaN
        return var
    elseif ndims==3
        tempvar=data.file[varname][:,:,:,itime]
        var=offset.+scale_factor.*tempvar
        var[tempvar.==dummy].=NaN
        return var
    else
        error("only 2D (one layer) and 3D (multiple layer) variables supported for now.")
    end
end

"""
   t_ref = get_reftime(Zarr_data)
Read the reference time from the attributes of time in the zarr file.
Times are described by a reference time (DateTime type) and the number of hours
relative to this t_ref. We use the midnight of the first time as t_ref.
"""
function get_reftime(data::ZarrData)
   time_relative=data.file.arrays["time"]
   units=time_relative.attrs["units"]
   first_time=CFTime.timedecode([time_relative[1]],units)
   # get year month day from first time
   yr=year(first_time[1])
   mo=month(first_time[1])
   dy=day(first_time[1])
   t_ref=DateTime(yr,mo,dy)
   # start of day of first time in dataset
   return t_ref
   # temp=split(units,"since")
   # temp2=split(temp[2],".")
   # println("ref date $(temp2[1])")
   # t0=DateTime(strip(temp2[1]),"yyyy-mm-dd HH:MM:SS")
   # dt_seconds=1.0
   # if startswith(temp[1],"seconds")
   #    dt_seconds=1.0
   # elseif startswith(temp[1],"minutes")
   #    dt_seconds=60.0
   # elseif startswith(temp[1],"hours")
   #    dt_seconds=3600.0
   # elseif startswith(temp[1],"days")
   #    dt_seconds=24.0*3600.0
   # else
   #    error("Invalid time-step unit in zarr-file.")
   # end
   # seconds_per_day=24.0*3600.0
   # if debuglevel>2
   #    println("$(dt_seconds), $(t0), $(seconds_per_day), $(time_relative[1])")
   # end
   # relative_days_in_seconds=div(dt_seconds*time_relative[1],seconds_per_day)*seconds_per_day
   # return t0 + Dates.Second(relative_days_in_seconds)
end

"""
   times=get_times(Zarr_data,Dates.DateTime(2019,1,1))

Get available times im netcdf map files as a range in seconds, eg 0.0:3600.0:7200.0
The reftime (second arg) is a Dates.DateTime can be used to change the reftime to something
convenient. With reftime=get_reftime(map) you can find the time of the first map, which is often
a nice default.
"""
function get_times(data::ZarrData,reftime::DateTime)
   time_relative=data.file.arrays["time"]
   units=time_relative.attrs["units"]
   times=CFTime.timedecode(time_relative[:],units)
   relative_times = [(t - reftime).value * 0.001 for t in times] # convert to seconds
   return relative_times
   # temp=split(units,"since")
   # temp2=split(temp[2],".")
   # println("ref date $(temp2[1])")
   # t0=DateTime(strip(temp2[1]),"yyyy-mm-dd HH:MM:SS")
   # dt_seconds=1.0
   # if startswith(temp[1],"seconds")
   #    dt_seconds=1.0
   # elseif startswith(temp[1],"minutes")
   #    dt_seconds=60.0
   # elseif startswith(temp[1],"hours")
   #    dt_seconds=3600.0
   # elseif startswith(temp[1],"days")
   #    dt_seconds=24.0*3600.0
   # else
   #    error("Invalid time-step unit in map-file.")
   # end
   # times=[]
   # for ti=1:length(time_relative)
   #    push!(times,(0.001*Dates.value(t0-reftime))+time_relative[ti]*dt_seconds)
   # end
   #return times
end

"""
   t = as_DateTime(data,reftime,relative_time)
Convert a reftime (DateTime) and relative_time (Float64 with seconds) relative to this.
Returns a DateTime type.
"""
function as_DateTime(data::ZarrData,reftime::DateTime,relative_time)
   return reftime+round(Int64,relative_time*1000)*Dates.Millisecond(1)
end

"""
p = initialize_interpolation(Zarr,"msl",t0)
Create an interpolation function p(x,y,z,t)
"""
function initialize_interpolation(data::ZarrData,varname::String,reftime::DateTime,dummy=0.0,cache_direction::Symbol=:forwards)
   times=get_times(data,reftime)
   values=data.file.arrays[varname] #TODO more checks
   missing_value=values.attrs["_FillValue"]
   scaling=values.attrs["scale_factor"]
   offset=values.attrs["add_offset"]
   xyt=CartesianXYTGrid(data.grid,times,values,varname,missing_value,scaling,offset,cache_direction)
   function f(x,y,z,t)
      value=interpolate(xyt,x,y,t,dummy)
      return value
   end
   return f
end
