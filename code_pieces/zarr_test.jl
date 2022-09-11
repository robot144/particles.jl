using Zarr
using NetCDF
using Dates
#include("src/dflow_his.jl")

# example 1 - one variable
z1 = zcreate(Int, 10000,10000,path = "example1.zarr",chunks=(1000, 1000))
z1[:] .= 42
z1[:,1] = 1:10000
z1[1,:] = 1:10000

z2 = zopen("example1.zarr")
z2[1:10,1:10]

#example 2
groupattrs = Dict("String attribute"=>"One", "Int attribute"=>5, "Float attribute"=>10.5)
g = zgroup("example2.zarr",attrs=groupattrs)

att1 = Dict("attribute"=>5)
v1 = zcreate(Int16, g, "var1_int",10,121,attrs=att1, chunks = (10,121))
v1[:,:].=21
att2 = Dict("attribute"=>15)
v2 = zcreate(Float32, g, "var2_float",9,10,11,attrs=att2, chunks = (9,10,11))
v2[:,:,:]=1.001*reshape(1:(9*10*11),(9,10,11))
att3 = Dict("attribute"=>15)
v3 = zcreate(Int16, g, "var3_int16",(5,6)...,attrs=att3, chunks = (10,10))
v3[:]=reshape(1:30,(5,6))

# example 3 - Strings
z3 = zcreate(Zarr.MaxLengthString{256,UInt8}, 3,path = "example3.zarr")
z3[:] = ["bla","blabla","blablabla"]

z3 = zopen("example3.zarr")
z3
# example 4 - remote dataset on nextcloud
z4=zopen("https://nx7384.your-storageshare.de/apps/sharingpath/wetwin/public/zunormm/ZUNO_his.zarr")
t4=z4["time"]


# inspecting his-file
s=NetCDF.open("test_data/locxxz_his.nc")

# inspecting map-file
m=NetCDF.open("test_data/locxxz_map.nc")
