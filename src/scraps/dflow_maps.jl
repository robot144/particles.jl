
using NetCDF
using Plots
Plots.default(:size,(1200,1200))

include("unstructured_grid.jl")
using Dates


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
   # mesh2d_face_nodes [4,:]
   #lat1=in1.vars["mesh2d_fourier001_min"]
   #xnodes1=in1.vars["NetNode_x"][:]
   for i=1:length(map)
      edges_temp=map[i].vars["NetElemNode"][:,:]
      xnodes_temp=map[i].vars["NetNode_x"][:]
      ynodes_temp=map[i].vars["NetNode_y"][:]
      @printf("- index computation\n")
      @time grid_temp=Grid(xnodes_temp,ynodes_temp,edges_temp,nmin,spherical)
      #dump(grid_temp)
      add_grid!(interp,grid_temp)
   end
   return interp
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
   t_start = (0.001*Dates.value(reftime-t0))+time_relative[1]*dt_seconds
   t_step  = (time_relative[2]-time_relative[1])*dt_seconds
   t_stop  = (0.001*Dates.value(reftime-t0))+time_relative[end]*dt_seconds
   times = t_start:t_step:t_stop
   return times
end

function as_DateTime(reftime::DateTime,relative_time)
   reftime+Float64(relative_time)*Dates.Second(1)
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


#
# main script
#

#dflow_map=load_nc_info("../test_data",r"rmm_toy_map.nc")
dflow_map=load_nc_info("../test_data",r"estuary_0000_map.nc")
interp=load_dflow_grid(dflow_map,50,false)
sealevel=load_nc_map_slice(dflow_map,"s1",10)

#global 
#xpoints=collect(range(-180.,stop=180.,length=300))
#ypoints=collect(range(-90.,stop=90.,length=300))
xpoints=collect(range(0.,stop=100000.,length=300))
ypoints=collect(range(0.,stop=500.,length=10))
sealevel_interp=interpolate(interp,xpoints,ypoints,sealevel)

using Plots
#contourf(xpoints,ypoints,sealevel_interp',clims=(-0.5,0.5)) #Plots
heatmap(xpoints,ypoints,sealevel_interp',clims=(-0.5,0.5)) #Plots
#heatmap(randn(100,100),size=(1000,1000),leg=false,border=:none,ticks=nothing)
# workaround padding for GR
#heatmap(randn(10,10),size=(500,500),leg=false,border=:none,ticks=nothing,left_margin=-4mm,bottom_margin=-4mm,top_margin=-2mm,right_margin=-2mm)
#pyplot()
#pcolormesh(xpoints,ypoints,lat_msl_interp',vmin=-1,vmax=1) #PyPlot
#savefig("fig_map_sealevel.png")

