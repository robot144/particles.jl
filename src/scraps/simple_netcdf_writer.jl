# Simple NetCDF writer for regular lat-lon grids
# 
using NetCDF

function nc_grid_write(filename,lon,lat,var_values,varname,dumval;title="UNKNOWN")
   lonatts = Dict("longname" => "Longitude",
                   "units"    => "degrees east")
   latatts = Dict("longname" => "Latitude",
                  "units"    => "degrees north")
   latdim = NcDim("lat", lat, lonatts)
   londim = NcDim("lon", lon, latatts)

   gatts=Dict("title"=>title)

   varatts=Dict("long_name" => varname,"missing_value"=>dumval, 
                "grid_mapping" =>"crs");
   #myvar = NcVar(varname, [londim, latdim], atts=varatts, t=Float32,compress=2,chunksize=(20,200))
   myvar = NcVar(varname, [londim, latdim], atts=varatts, t=Float32)

   crsatts=Dict("grid_mapping_name"=>"latitude_longitude");
   crsvar = NcVar("crs", NcDim[], atts=crsatts, t=Int32)

   nc = NetCDF.create(filename, NcVar[myvar,crsvar],gatts=gatts,mode=NC_NETCDF4)

   NetCDF.putvar(nc, varname, var_values)
   crsval= Array{Int32,0}(undef);crsval[1]=4326
   NetCDF.putvar(nc, "crs", crsval)

   NetCDF.close(nc)
end

#
# testing
#
function test1()
   cfn2 = "ff2.nc"
   isfile(cfn2) && rm(cfn2)
   lat = collect(-89.0:1.0:89.0)
   lon = collect(-180.0:1.0:180.0)
   vals=randn(Float32,length(lon),length(lat))
   dumval=Float32(999.0)
   vals[10:20,10:40].=dumval
   nc_grid_write(cfn2,lon,lat,vals,"random",dumval,title="Created by me.")
end

#test1()
