# csv file with location over time of one of a set of buoys installed on sea-ice in the Arctic. The main 
# transport is southward just east of Svalbard.
# Bouy 16 : 2014 March 15 - May 15

using DataFrames
using CSV
using Plots
using Dates
using Particles

# read buoy track
b16=CSV.read("CliSAP_Boje_16_modified.csv")
# describe(b16)
b16_n=length(b16.Lon)
b16_times=[]
for ti=1:b16_n
   push!(b16_times,DateTime(b16[ti,1],"dd-mm-yyyy HH:MM"))
end
b16.Times=b16_times

# read GTSM current data
dflow_map=load_nc_info(".",r"gtsm_fine_...._map.nc")
t0=get_reftime(dflow_map) # March 1 2014
#tmaps=get_times(dflow_map,t0) #all times in file
# create interpolation functions u1(x,y,z,t) and v1(x,y,z,t) with t in seconds relative to t0
u1,v1=initialize_interpolation(dflow_map,interp,t0);
# some tests
#interp=load_dflow_grid(dflow_map,50,false)
#ind=find_index(interp,-2.0,50.0)
#ff=u1(22.0,75.0,0.0,) #lon,lat,dummy,t with t in seconds relative to t0
# create timeseries at one point
lon=22.0
lat=75.0
ts=(14.0*24.0*3600.0):1800.0:(21.0*24.0*3600.0) #times in seconds since t0
tt=t0+ts.*Second(1) #convert to DateTime
uu=zeros(length(tt))
vv=zeros(length(tt))
for ti=1:length(tt)
  uu[ti]=u1(lon,lat,0.0,ts[ti])
  vv[ti]=v1(lon,lat,0.0,ts[ti])
end
plot(tt,[uu,vv],label=["u east","v north"])
title!("Tidal currents from GTSM model at lon=22E lat=75N")
savefig("timseries_uv_22e_75n.png")
# first interpolate to track positions and times for buoy 16
utrack16=zeros(length(b16.Times))
vtrack16=zeros(length(b16.Times))
ttrack16=zeros(length(b16.Times))
println("Please, wait. It takes a few minutes to move through the 16Gb input file")
for ti=1:(length(ttrack16)-96) #!!! the last day is missing in the GTSM data
   time=(b16.Times[ti]-t0).value/1000 #convert to seconds relative to t0
   #println("time=$(time) $(ti)")
   ttrack16[ti]=time
   utrack16[ti]=u1(b16.Lon[ti],b16.Lat[ti],0.0,ttrack16[ti])
   vtrack16[ti]=v1(b16.Lon[ti],b16.Lat[ti],0.0,ttrack16[ti])
end

#
# plot bouy 16 on a background
#
width=1000
height=1000
Plots.default(:size,[width,height])
gebco_server=WmsServer("gebco")
bbox=[0.0,70.0,40.0,85.0]
img=get_map(gebco_server,bbox,width,height)
plot_image(img,bbox)
plot!(b16.Lon,b16.Lat,label=["Buoy 16"])
savefig("track_buoy16.png")

#
# plot GTSM tidal currents along track
#
plot(b16.Times,[utrack16,vtrack16],label=["u east","v north"])
ylabel!("Tidal current [m/s]")
title!("Tidal currents from GTSM along track of buoy 16")
savefig("tidal_currents_buoy16.png")

