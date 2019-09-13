# Interact with dflow files
#

using NetCDF
using Plots
using TimeSeries
using Dates

debug=0

"""
   true_or_false = his_hasvar("mynetcdf.nc","waterlevel")

Check for presence of a variable in a dflow netcdf history file.
"""
function his_hasvar(fromfile::String,varname::String)
   info=ncinfo(fromfile)
   return haskey(info.vars,varname)
end
function his_hasvar(fromfile::NcFile,varname::String)
   return haskey(fromfile.vars,varname)
end

"""
   times = his_get_times("mynetcdf.nc")
"""
function his_get_times(fromfile::String)
   info=ncinfo(fromfile)
   return his_get_times(info)
end
function his_get_times(fromfile::NcFile)
   time_relative=fromfile.vars["time"]
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
      error("Invalid time-step unit in his-file.")
   end
   #times=(t0+time_relative[1]*dt_seconds*Second(1)):((time_relative[2]-time_relative[1])*dt_seconds*Second(1)):(t0+time_relative[end]*dt_seconds*Second(1))
   times=[ (t0+time_relative[i]*dt_seconds*Second(1)) for i=1:length(time_relative) ]
   return times
end

"""
   stations = his_get_stations("mynetcdf.nc")
"""
function his_get_stations(fromfile::String)
   info=ncinfo(fromfile)
   return his_get_stations(info)
end
function his_get_stations(fromfile::NcFile)
   stations=nc_char2string(fromfile.vars["station_name"][:,:])
   return stations
end


# initial test
function test1()
   hisfile="../test_data/estuary_his.nc"
   @assert his_hasvar(hisfile,"waterlevel")
   @assert his_hasvar(hisfile,"notavar")==false
   t=his_get_times(hisfile)
   @assert t[1]==Dates.DateTime(1991,1,1,0,0,0)
   @assert t[34561]==Dates.DateTime(1991,8,29,0,0,0)
   
   stations=his_get_stations(hisfile)
   @assert stations[2]=="station02"
   @assert findfirst(x->isequal(x,"station03"),stations)==3
   @assert findfirst(x->isequal(x,"station04"),stations)==nothing
end
test1()
