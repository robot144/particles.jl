# simple simulation to study the behaviour of floating ice 
# under the influence of winds
using Plots
seconds= 1.0 #si unit
minutes= 60.0*seconds
hours  = 60.0*minutes  
days   = 24.0*hours
rad2deg= 180.0/pi
deg2rad= pi/180.0

#simulation parameters
tstart = 0.0 #arbitrary initial time
tstep  = 10*minutes
tstop  = 10.0*days
tstep  = 10.0*minutes
lat    = 76.0 
omega  = (2.0*pi)/(24.0*hours)
f      = 2.0*omega*sin(lat*deg2rad)
r_ice_water = 5.5e-3
r_air_ice = 1.2e-3
rho_ice = 920.0
rho_air = 1.2
rho_water = 1020.0
u_wind_max = 15.0 # m/s in x-direction
u_wind_min = 1.0
h = 3.0 #ice thickness

function u_wind(t)
   result=u_wind_min
   if (t>(6.25*days))
      result=u_wind_min
   elseif (t>(6*days))
      result = (u_wind_max-u_wind_min)*(1.0-(t-(6*days))/(0.25*days))+u_wind_min
   elseif (t>(5*days))
      result=u_wind_max
   elseif (t>(4*days))
      result = (u_wind_max-u_wind_min)*(t-(4*days))/(1*days)+u_wind_min
   end
   return result
end
a_tide = 0.3 #tidal amplitude [m]
g=9.81
D=500.0 #depth [m]
omega_m2=(2*pi)/(12.4*hours)
function u_tide(t)
   return a_tide*sqrt(g/D)*sin(omega_m2*t)
end

# model dx/dt = model(x,t)
# a = F/m Newton budget per square meter of surface area
function model(x,t)
   xc = x[1]
   yc = x[2]
   u  = x[3]
   v  = x[4]
   dx=zeros(4)
   dx[1] = d_xc = u
   dx[2] = d_yc = v
   u_mag = sqrt((u-u_tide(t))^2+v^2)
   dx[3] = d_u  = (r_ice_water*rho_water*(u_tide(t)-u)*u_mag + 
		   r_air_ice*rho_air*u_wind(t)^2)/(rho_ice*h) +f*v
   dx[4] = d_v  = (r_ice_water*rho_water*(-v)*u_mag)/(rho_ice*h) -f*u
   return dx
end

#run simulation
times = tstart:tstep:tstop
x=[0.0, 0.0, 0.0, 0.0] #initial values
results=zeros(length(x),length(times))
ti=1
for t=times
   global x
   global ti
   #println("t=$(t) x=$(x)")
   results[:,ti]=x[:]
   xn = x + tstep*model(x,t) 
   x  = xn
   ti+=1
end

#plot results
px=plot(times/hours,results[1,:],ylabel="x")
py=plot(times/hours,results[2,:],ylabel="y")
pu=plot(times/hours,results[3,:],ylabel="u")
pv=plot(times/hours,results[4,:],ylabel="v")
plot(px,py,pu,pv,layout=(4,1),size=(1200,800),legend=false)
#savefig("fig_series_ice.png")

#plot(results[1,:],results[2,:],xlabel="x",ylabel="y",size=(1200,1200))
#savefig("fig_track_ice.png")

