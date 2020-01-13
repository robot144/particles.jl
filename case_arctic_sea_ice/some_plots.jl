#
# create some plots
#
using Particles
using NetCDF
using Plots

# CMEMS data
cmems=CmemsData(".","cmems_201403_05.nc")
x=cmems.grid.xnodes
y=cmems.grid.ynodes
x=range(x[1],x[end],length=length(x)) #remove slight differences in grid size
y=range(y[1],y[end],length=length(y))
t0=get_reftime(cmems)
times=get_times(cmems,t0)
ti=1
uice=load_map_slice(cmems,"usi",ti) #sea ice velocity eastward
vice=load_map_slice(cmems,"vsi",ti) #sea ice velocity northward
cice=load_map_slice(cmems,"siconc",ti) #sea ice concentration
hice=load_map_slice(cmems,"sithick",ti) #sea ice thickness

uice_mag=sqrt.(uice.^2 + vice.^2)
isice=(hice.>0.2)
landice=isice .& (uice_mag.<0.005)
temp=zeros(size(landice))
temp[.!landice].=NaN

#background
width=1000
height=1000
Plots.default(:size,[width,height])
gebco_server=WmsServer("gebco")
bbox=[0.0,70.0,40.0,85.0]
img=get_map(gebco_server,bbox,width,height)

#estimate of landfast ice
plot_image(img,bbox)
heatmap!(x,y,temp')

