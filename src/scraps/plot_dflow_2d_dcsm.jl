#
# read output file of run and create some plots
#
using Plots
width=2000
height=1000
Plots.default(:size,[width,height])

using NetCDF
include("wms_client.jl") 

#ncfilename="output_dflow_2d_dcsm.nc"
ncfilename="output_dflow_2d_dcsm2019.nc"
#fig_prefix="fig_dflow_2dcsm"
fig_prefix="fig_dflow_2dcsm2019"

ncfile=ncinfo(ncfilename)
lon=ncfile["lon"][:,:]
lat=ncfile["lat"][:,:]

wms_server=WmsServer("emodnet-bathymetry")
#bbox=[-10.0,48.0,10.0,60.0]
bbox=[4.5,53.0,5.5,54.0]
img=get_map(wms_server,bbox,width,height)


# initial position of particles
plot_image(img,bbox)
scatter!(lon[:,1],lat[:,1])
savefig("$(fig_prefix)_initial.png")

# final position of particles
plot_image(img,bbox)
scatter!(lon[:,end],lat[:,end])
savefig("$(fig_prefix)_final.png")

# final position of particles
plot_image(img,bbox)
plot!(lon[:,:]',lat[:,:]')
savefig("$(fig_prefix)_tracks.png")
