# Flow field

using Plots
using StaticArrays
#using Profile

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
      println("t=$(t)\n")
      for i=1:n
         s[:]=p[:,i]
         f!(ds,s,t,i,d)
         s+=ds*dt #Euler forward
         p[:,i]=s[:]
      end
      t+=dt
   end
   return t
end

function run_simulation(d)
   Plots.default(:size,d["plot_maps_size"])

   if d["plot_maps"]
      global fig1=d["plot_maps_background"](d)
      #scatter!(fig1,p[1,:],p[2,:],legend=false)
      d["plot_maps_func"](fig1,d,p)
      # TODO possibly label=[string(t)])
      gui(fig1)
   end
   
   t=d["tstart"]
   tend=d["tend"]
   target_times=d["plot_map_times"]
   if((length(target_times)==0) || (target_times[end]<tend))
	push!(target_times,tend)
   end
   for t_stop=target_times
      simulate!(p,t,t_stop,d)
      t=t_stop
      if d["plot_maps"]
         #scatter!(fig1,p[1,:],p[2,:])
         d["plot_maps_func"](fig1,d,p)
         #sleep(1)
         gui(fig1)
      end
   end
   
   #wait for user 
   if !isinteractive() #wait for user to kill final plot
      println("Type [enter] to finish script")
      readline()
   end
end

function index(var,vars)
   return indexin([var],vars)[1]
end

function plot_maps_xy(fig,d,p)
   x_index=index("x",d["variables"])
   y_index=index("y",d["variables"])
   scatter!(fig,p[x_index,:],p[y_index,:],markercolor=[:black],legend=false)
end

function plot_maps_xz(fig,d,p)
   x_index=index("x",d["variables"])
   z_index=index("z",d["variables"])
   scatter!(fig,p[x_index,:],p[z_index,:],legend=false)
end

function default_userdata()
   d=Dict()
   d["dt"]=0.01 #TODO a fixed timestep will not work in general
   d["tstart"]=0.0
   d["tend"]=1.0
   d["plot_maps"]=true
   d["plot_maps_size"]=(1000,1000)
   d["plot_maps_times"]=[]
   d["plot_maps_func"]=plot_maps_xy
   return d
end

