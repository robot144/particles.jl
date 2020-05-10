
# make a plot of the era5 data

using Plots
using Dates
using Particles

# background
width=1000
height=1000
Plots.default(:size,[width,height])
gebco_server=WmsServer("gebco") #gebco or emodnet-bathymetry or open-streetmap
plotbox=[5.0,53.0,11.0,56.0] #area to plot min(Lon), min(Lat), max(Lon), max(Lat)

img=get_map(gebco_server,plotbox,width,height)

plot_image(img,plotbox)

#era5 data
era=EraData(".","era5_wind_201703_04.nc")
x=era.grid.xnodes
y=era.grid.ynodes
t0=get_reftime(era)
times=get_times(era,t0)

itime=1
t_rel=times[itime]
t=as_DateTime(era,t0,t_rel)
println("time = $(t)")

hs=load_map_slice(era,"swh",itime)
#heatmap!(x,y,hs')
# !!! workaround is flip for plotting
heatmap!(x,y[end:-1:1],hs[:,end:-1:1]')
title!("Hs at $(t)")

#!!!! Issue y is in DECREASING order !!!


# function initialize_interpolation(data::EraData,varname::String,reftime::DateTime,dummy=0.0)

