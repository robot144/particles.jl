
# make a plot of the swan data

using Plots
using Dates
using Particles

# background
width=1000
height=1000
Plots.default(:size,[width,height])
gebco_server=WmsServer("gebco") #gebco or emodnet-bathymetry or open-streetmap
plotbox=[3.0,50.0,13.0,56.0] #area to plot min(Lon), min(Lat), max(Lon), max(Lat)

img=get_map(gebco_server,plotbox,width,height)


#swan data from matroos
swan=MatroosData(".","swan_20171007_28.nc")
x=swan.grid.xnodes
y=swan.grid.ynodes
t0=get_reftime(swan)
times=get_times(swan,t0)

itime=57
itime=80
t_rel=times[itime]
t=as_DateTime(swan,t0,t_rel)
println("time = $(t)")

hs=load_map_slice(swan,"wave_height_hm0",itime)

plot_image(img,plotbox)
heatmap!(x,y,hs')
title!("Hs at $(t)")


# function initialize_interpolation(data::EraData,varname::String,reftime::DateTime,dummy=0.0)

