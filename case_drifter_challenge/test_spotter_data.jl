

using DataFrames
using CSV
using Plots
using Particles
using Dates

spotterfile="data/hindcast/challenge_30-days_sofar_20210530_atlantic_sample.csv"

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

# collect positions at initial time
# common reference is 2021-05-30T00:15:00 but the first time per spotter varies a bit
t0=DateTime("2021-05-30T00:15:00")
lon=zeros(ndrifters)
lat=zeros(ndrifters)
id=collect(0:(ndrifters-1))
for i=1:ndrifters
   println("i=$(i)")
   dfi = df[df.spotId.==(i-1),[:time,:lon,:lat]]
   w=(t0-dfi.time[1])/(dfi.time[2]-dfi.time[1])
   lon[i]=(1-w)*dfi.lon[1]+w*dfi.lon[2]
   lat[i]=(1-w)*dfi.lat[1]+w*dfi.lat[2]
   if abs(w)>2
      error("drifter $(i) does not start close to t0. $(w)")
   end
end


width=1000
height=600
gebco_server=WmsServer("gebco")
boundbox=[-90.0,20.0,-10.0,50.0]
img=get_map(gebco_server,boundbox,width,height)
plot_image(img,boundbox)
plot!(lon,lat,seriestype = :scatter, label = "initial position")
savefig("initial_spotter_positions.png")
