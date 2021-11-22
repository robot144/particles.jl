

using DataFrames
using CSV
using Plots
using Particles
using Dates
using NetCDF
import Interpolations

R = 6371.0e3               # Mean radius of earth from wikipedia

#
# hindcast
#

spotterfile="data/hindcast/challenge_30-days_sofar_20210530_atlantic_sample.csv"
particlefile="netcdf_output/darpa_hcast.nc"
#runid="leeway0.0144npart1k1.0"
#runid="leeway0.0072npart1k1.0"
#runid="leeway0.0288npart1k1.0"
#runid="leeway0.0100npart1k1.0"
runid="leeway0.0144npart10k500.0"

# load tracks and look for starting points
df = DataFrame(CSV.File(spotterfile))
#simplify some column names
df.lon=select(df,"longitude [decimal degrees]")[:,1] .- 360.0
df.lat=select(df,"latitude [decimal degrees]")[:,1]
df.timeStrings=select(df,"timestamp [utc]")[:,1]
n=size(df,1)
times=zeros(DateTime,n)
for ti=1:n
   times[ti]=DateTime(df[ti,:timeStrings],"yyyy-mm-dd HH:MM:SS+00:00")
end
df.time=times
ndrifters=maximum(df.spotId)+1 # starting with spotId 0

# create simplified data structure with the relevant measurements
measurements=[]
t0=DateTime("2021-05-30T00:15:00")
lon=zeros(ndrifters)
lat=zeros(ndrifters)
id=collect(0:(ndrifters-1))
for i=1:ndrifters
   println("i=$(i)")
   dfi = df[df.spotId.==(i-1),[:time,:lon,:lat]]
   push!(measurements,dfi)
   w=(t0-dfi.time[1])/(dfi.time[2]-dfi.time[1])
   lon[i]=(1-w)*dfi.lon[1]+w*dfi.lon[2]
   lat[i]=(1-w)*dfi.lat[1]+w*dfi.lat[2]
   if abs(w)>2
      error("drifter $(i) does not start close to t0. $(w)")
   end
end

# read model output
output=NetCDF.open(particlefile)
p_lon=output.vars["lon"][:,:,1]
p_lat=output.vars["lat"][:,:,1]
n=size(p_lon,2) # total nummber of particles
npart=div(n,ndrifters)
tref=DateTime("2021-05-30T00:15:00") #TODO read from netcdf
p_time=output.vars["time"][:]
istart=2
# drifter 12 has a duplicate time at nr approx 492 -> removed line manually from input
# same drifter 59 around 72
#iend=550
# interpolate in time and compute errors
iend=490
t_meas  =p_time[istart:iend]
ntime=length(t_meas)
rmse=zeros(ndrifters,ntime)
p_lon_avg=zeros(ntime,ndrifters)
p_lat_avg=zeros(ntime,ndrifters)
for i=1:ndrifters
   println("drifter $(i)")
   t_tmp=div.(measurements[i].time-tref,Dates.Millisecond(1000))
   lon_tmp=measurements[i].lon
   lat_tmp=measurements[i].lat
   flon=Interpolations.LinearInterpolation(t_tmp,lon_tmp)
   flat=Interpolations.LinearInterpolation(t_tmp,lat_tmp)
   lon_meas=flon.(t_meas) #index 1 outside interval for measurements
   lat_meas=flat.(t_meas)
   imin=npart*(i-1)+1
   imax=npart*i
   lon_mod_all=p_lon[istart:iend,imin:imax]
   lat_mod_all=p_lat[istart:iend,imin:imax]
   lon_mod=sum(lon_mod_all,dims=2)/size(lon_mod_all,2)
   lat_mod=sum(lat_mod_all,dims=2)/size(lat_mod_all,2)
   p_lon_avg[:,i].=lon_mod[:]
   p_lat_avg[:,i].=lat_mod[:]
   dx=1e-3*R*deg2rad.(lon_mod-lon_meas).*cosd.(lat_meas)
   dy=1e-3*R*deg2rad.(lat_mod-lat_meas)
   dist=sqrt.(dx.^2 .+ dy.^2)
   rmse[i,:].=dist[:]
end

out=open("rmse_$(runid).json","w")
JSON.print(out,Dict("t_meas"=>t_meas,"rmse"=>rmse))
close(out)


#rmse[77,:].=0.0 # has NaN values. It is probably on land in the model
rmse[isnan.(rmse)].=0.0 #affects stats but is needed
rmse_mean=sum(rmse,dims=1)/size(rmse,1)

# plot stats
plot(t_meas/3600.0/24.0,rmse',color=:grey,label="")
plot!(t_meas/3600.0/24.0,rmse_mean',color=:red,lw=3,label="RMSE")
xlabel!("lead-time [days]")
ylabel!("RMSE distance [km]")
title!("runid=$(runid)")
ylims!((0.0,500.0))
savefig("distance_rmse_$(runid).png")

# plot some tracks
width=2000
height=1000
gebco_server=WmsServer("gebco")
boundbox=[-90.0,20.0,-10.0,50.0]
img=get_map(gebco_server,boundbox,width,height)
plot_image(img,boundbox)
plot!(lon,lat,seriestype = :scatter, label = "initial position")
for i=1:length(measurements)
   plot!(measurements[i].lon,measurements[i].lat,color=:red,label="")
end
title!("Measured tracks")
xlabel!("longitude")
ylabel!("latitude")
savefig("measured_spotter_positions.png")

plot_image(img,boundbox)
plot!(p_lon,p_lat,color=:lightblue,label="")
plot!(lon,lat,seriestype = :scatter, lw=2,label = "initial position")
for i=1:length(measurements)
   plot!(measurements[i].lon,measurements[i].lat,color=:red,label="")
end
plot!(p_lon_avg,p_lat_avg,color=:blue,lw=2,label="")
title!("Modeled and measured tracks")
xlabel!("longitude")
ylabel!("latitude")
savefig("modelled_spotter_positions_$(runid).png")

