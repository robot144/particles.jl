# csv file with location over time of one of a set of buoys installed on sea-ice in the Arctic. The main 
# transport is southward just east of Svalbard.
# Bouy 16 : 2014 March 15 - May 15

using DataFrames
using CSV
using Plots
using Dates
using Particles

# plot a background
width=1000
height=1000
Plots.default(:size,[width,height])
gebco_server=WmsServer("gebco")
bbox=[0.0,70.0,40.0,85.0]
img=get_map(gebco_server,bbox,width,height)
plot_image(img,bbox)

# add buoy track
b16=CSV.read("CliSAP_Boje_16_modified.csv")
lon=b16[:Lon]
lat=b16[:Lat]

plot!(lon,lat,label=["Buoy 16"])
#savefig("track_buoy16.png")
