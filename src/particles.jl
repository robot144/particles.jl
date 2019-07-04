# Shared routines for particle modelling. 
# No routines for specific data sources or particle model-equations should be included,
# just the generic stuff goes here.

using Plots
using StaticArrays
using Dates

"""
   t_next=simulate!(p,t_now,t_stop,d)

Compute multiple timesteps until t>=t_stop for all particles in matrix p.
Possibly using parameters (eg dt) from d of type userdata.
"""
function simulate!(p,t,t_stop,d)
   dt=d["dt"]
   f=d["f"]
   variables=d["variables"]
   (m,n)=size(p) # no variables x no particles
   ds=@MVector zeros(length(variables))
   s=@MVector zeros(length(variables))
   while(t<t_stop)
      #println("... t=$(t) < $(t_stop)")
      for i=1:n
         s[:]=p[:,i]
         f!(ds,s,t,i,d)
         s+=ds*dt #Euler forward
         #println("   $(i) $(s)")
         p[:,i]=s[:]
      end
      t+=dt
   end
   return t
end

"""
   run_simulation(d)
   Main simulation routine. The configuration is contained in Dict d, 
   which may contain call-backs for plotting etc.
"""
function run_simulation(d)
   vars=d["variables"]
   npart=d["nparticles"]
   nvars=length(vars)
   Plots.default(:size,d["plot_maps_size"])

   if d["plot_maps"]
      Plots.default(:size,d["plot_maps_size"])
      global fig1=d["plot_maps_background"](d)
      #scatter!(fig1,p[1,:],p[2,:],legend=false)
      d["plot_maps_func"](fig1,d,p)
      # TODO possibly label=[string(t)])
      gui(fig1)
   end
   if d["write_maps"]
      (nc_out,ncvars)=initialize_netcdf_output(d)
   end
   
   t=d["tstart"]
   tend=d["tend"]
   plot_maps_times=d["plot_maps_times"]
   write_maps_times=d["write_maps_times"]
   #break up the simulation in chunks that can run without producing output
   target_times=sort(union(plot_maps_times,write_maps_times))
   if((length(target_times)==0) || (target_times[end]<tend))
	push!(target_times,tend)
   end
   if target_times[end]>tend
      temp=ssort(union(target_times,tend))
      i_last=findlast(x->x<=tend,temp)
      target_times=temp[1:i_last]
   end
   # simulate in pieces until next output-action
   for t_stop=target_times
      t=simulate!(p,t,t_stop,d)
      println("t=$(t)")
      if (d["plot_maps"]) && (t_stop in plot_maps_times)
         #scatter!(fig1,p[1,:],p[2,:])
         d["plot_maps_func"](fig1,d,p)
         #sleep(1)
         gui(fig1)
      end
      if (d["write_maps"]) && (t_stop in write_maps_times)
         timei=findfirst(x->x==t_stop,write_maps_times)
         for vari=1:nvars
            varname=vars[vari]
            start=[1,timei] #part x time
            count=[npart,1]
            NetCDF.putvar(ncvars[vari], p[vari,:]; start=start, count=count)
         end
      end
   end

   NetCDF.close(nc_out)
   
   #wait for user 
   if !isinteractive() #wait for user to kill final plot
      println("Type [enter] to finish script")
      readline()
   end
end

"""
   i1 = index(2,[1,2,3,4])
   i2 = index("bob",["alex","bob","charlie"])
   Find first occurrence of a variable in an array. 
   Returns nothing if the value is not found.
"""
function index(var,vars)
   return indexin([var],vars)[1]
end

"""
   f=plot_maps_xy(fig,d,p)
   Plot particles as dots in xy-plane.
"""
function plot_maps_xy(fig,d,p)
   x_index=index("x",d["variables"])
   y_index=index("y",d["variables"])
   scatter!(fig,p[x_index,:],p[y_index,:],markercolor=[:black],legend=false)
end

"""
   f=plot_maps_xz(fig,d,p)
   Plot particles as dots in xy-plane.
"""
function plot_maps_xz(fig,d,p)
   x_index=index("x",d["variables"])
   z_index=index("z",d["variables"])
   scatter!(fig,p[x_index,:],p[z_index,:],legend=false)
end

function initialize_netcdf_output(d)
   # file 
   filedir=d["write_maps_dir"]
   filename=d["write_maps_filename"]
   fullfile=joinpath(filedir,filename)
   if !isdir(filedir)
      println("Directory for output does not exist $(filedir). Creating it.")
      mkpath(filedir)
   end
   if isfile(fullfile)
      println("Output file exists. Removing file $(fullfile)")
      rm(fullfile)
   end
   # time
   t0=d["reftime"]
   time_atts = Dict("units" => "seconds since $(Dates.format(t0,"yyyy-mm-dd HH:MMM:SS"))",
                    "standard_name" => "time", "long_name" => "time",
		    "comment" => "unspecified time zone", "calendar" => "gregorian" )
   map_times=collect(d["write_maps_times"])
   time_dim = NcDim("time", map_times, time_atts)
   # particle
   npart=d["nparticles"]
   vars=d["variables"]
   nvars=length(vars)
   part_atts = Dict("long_name" => "particle id")
   part_dim = NcDim("particles", collect(1:1:npart), part_atts)

   #global attributes
   gatts=Dict("title"=>"Output of particle model")

   myvars=[]
   dumval=d["dumval"]
   for vari=1:nvars
      varname=vars[vari]
      varatts=Dict("long_name" => varname,"missing_value"=>dumval) 
      if varname=="x"
          varatts["long_name"] = "x-coordinate"
          varatts["units"]     = "m"
      elseif varname=="y"
          varatts["long_name"] = "y-coordinate"
          varatts["units"]     = "m"
      elseif varname=="lon"
          varatts["long_name"] = "Longitude"
          varatts["units"]     = "degrees east"
      elseif varname=="lat"
          varatts["long_name"] = "Latitude"
          varatts["units"]     = "degrees north"
      end
      myvar = NcVar(varname, [part_dim, time_dim], atts=varatts, t=Float32)

      push!(myvars,myvar)
   end

   nc = NetCDF.create(fullfile, NcVar[myvars...],gatts=gatts,mode=NC_NETCDF4)

   p=d["particles"] # var x part
   for vari=1:nvars
      varname=vars[vari]
      start=[1,1] #part x time
      count=[npart,1]
      NetCDF.putvar(myvars[vari], p[vari,:]; start=start, count=count)
   end

   return (nc,myvars)
   #NetCDF.close(nc)
end

"""
   d = default_userdata()
   initialize Dict with some default configuration data
"""
function default_userdata()
   d=Dict()
   # general
   d["dt"]=0.01 #TODO a fixed timestep will not work in general
   d["tstart"]=0.0
   d["tend"]=1.0
   d["reftime"]=DateTime(2000,1,1) # Jan 1st 2000
   d["coordinates"]="projected" #projected or spherical 
   d["nparticles"]=10 #number of particles
   d["variables"]=["x","y"] #recognized are x,y,lat,lon other variables are written with partial meta-data
   d["dumval"]=9999.0
   #plotting to screen
   d["plot_maps"]=false
   d["plot_maps_size"]=(1200,1000)
   d["plot_maps_times"]=[]
   d["plot_maps_func"]=plot_maps_xy
   #results that are kept in memmory
   d["keep_particles"]=false
   d["keep_particle_times"]=[]
   #results written to netcdf file
   d["write_maps"]=false
   d["write_maps_times"]=[]
   d["write_maps_dir"]="."
   d["write_maps_filename"]="output.nc"
   return d
end

