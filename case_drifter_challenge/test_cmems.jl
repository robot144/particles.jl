# read some cmems data for testing


#
# explore with NetCDF module
#
using NetCDF
using Plots
u=NetCDF.open("data/hindcast/u0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc")
u_scale=u.vars["uo"].atts["scale_factor"]
u_dummy=u.vars["uo"].atts["_FillValue"]
v=NetCDF.open("data/hindcast/v0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc")
v_scale=v.vars["vo"].atts["scale_factor"]
v_dummy=v.vars["vo"].atts["_FillValue"]
lon=u.var["longitude"]
lat=u.var["latitude"]
itime=1 # time index
uint=u.vars["uo"][:,:,1,itime]
ui=u_scale*uint
ui[uint.==u_dummy].=NaN
vint=v.vars["vo"][:,:,1,itime]
vi=v_scale*uint
vi[vint.==v_dummy].=NaN
heatmap(lon,lat,sqrt.(ui.^2+vi.^2)')


#
# read with cmems_grid.jl
#
using Particles

cmems_u = CmemsData("data/hindcast","u0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc")
cmems_v = CmemsData("data/hindcast","v0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc")

u1=load_map_slice(cmems_u,"uo",1) # u-velocity

t0=get_reftime(cmems_u)

times=get_times(cmems_u,t0)

t2=as_DateTime(cmems_u,t0,times[2])

u=initialize_interpolation(cmems_u,"uo",t0,NaN) #water velocity x-dir
v=initialize_interpolation(cmems_u,"vo",t0,NaN) #water velocity y-dir
u1=u(-73.0,35.0,0.0,1800.0)
v1=v(-73.0,35.0,0.0,1800.0)

