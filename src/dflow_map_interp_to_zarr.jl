#
# Convert delft3d time-series to zarr format
# This add selection, compression and chuning into separate files per chunk
#
using TOML
using NetCDF
using Zarr
using Dates
if !@isdefined(Grid) #TODO: This is awkward. Should move code to a package 
    include(joinpath(@__DIR__,"unstructured_grid.jl"))
    include(joinpath(@__DIR__,"dflow.jl"))
end

debuglevel=1

#
# defaults
#
# these variables are added to the configuration by default
try_vars = ["waterlevel","x_velocity","y_velocity","z_velocity","salinity","temperature","waterdepth",
    "eddy_visc_z", "z_center_3d","z_iface_3d"] 
# default settings per variable
defaults = Dict(
    "waterlevel" => Dict(
        "name" => "waterlevel",
        "scale_factor" => 0.001, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "x_velocity" => Dict(
        "name" => "x_velocity",
        "scale_factor" => 0.001, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "y_velocity" => Dict(
        "name" => "y_velocity",
        "scale_factor" => 0.001, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "z_velocity" => Dict(
        "name" => "z_velocity",
        "scale_factor" => 1.0E-6, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "z_velocity_center" => Dict(
        "name" => "z_velocity_center",
        "scale_factor" => 1.0E-6, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "salinity" => Dict(
        "name" => "salinity",
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "temperature" => Dict(
        "name" => "temperature",
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "eddy_visc_z" => Dict(
        "name" => "eddy_visc_z",
        "scale_factor" => 1.0E-6, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "waterdepth" => Dict(
        "name" => "waterdepth",
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "time" => Dict(
        "scale_factor" => 1.0, 
        "add_offset" => 0.0,
        "data_type" => "Float64",
        "_FillValue" => -9999.0),
    "x_center" => Dict(
        "scale_factor" => 1.0, 
        "add_offset" => 0.0,
        "data_type" => "Float64",
        "_FillValue" => -9999.0),
    "y_center" => Dict(
        "scale_factor" => 1.0, 
        "add_offset" => 0.0,
        "data_type" => "Float64",
        "_FillValue" => -9999.0),
    "z_center" => Dict(
        "scale_factor" => 1.0, 
        "add_offset" => 0.0,
        "data_type" => "Float64",
        "_FillValue" => -9999.0),
    "z_center_3d" => Dict( #z at center (depth limited to around 3000m)
        "scale_factor" => 0.1, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => -9999.0),
    "z_iface_3d" => Dict( #z at interface (depth limited to around 3000m)
        "scale_factor" => 0.1, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => -9999.0)
         )

# variables appear under different names in the delft3d-fm output files. Here we list the options
aliases=Dict{String,Vector{String}}(
    "waterlevel"  => ["s1", "mesh2d_s1"],
    "x_velocity"  => ["ucx", "mesh2d_ucx"],
    "y_velocity"  => ["ucy", "mesh2d_ucy"],
    "z_velocity"  => ["ww1", "mesh2d_ww1"], #vertical velocity at interfaces
    "z_velocity_center"  => ["ucz", "mesh2d_ucz"], #vertical velocity at layer centers
    "salinity"    => ["sa1","mesh2d_sa1"], #sa1 not sal (one not L)?
    "waterdepth"  => ["waterdepth","mesh2d_waterdepth"],
    "temperature" => ["tem1","mesh2d_tem1"], 
    "eddy_visc_z" => ["vicwwu","mesh2d_vicwwu"], 
    "x_center"    => ["FlowElem_xcc","mesh2d_face_x"],
    "y_center"    => ["FlowElem_ycc","mesh2d_face_y"],
    "z_center"    => ["mesh2d_layer_z","LayCoord_cc"], #1d
    "z_center_3d" => ["mesh2d_flowelem_zcc"], #3d
    "z_iface_3d"  => ["mesh2d_flowelem_zw"], #3d
    "x_node"      => ["mesh2d_node_x","NetNode_x"],
    "y_node"     => ["mesh2d_node_y","NetNode_y"],
    "time"       => ["time"]
)
# valriable names that have a vertical position at the interface instead of the center. The center is the defauls
znames_iface=["z_iface_3d","z_velocity", "eddy_visc_z"]

chunk_target_size=1000000
#
# supporting functions
#
typeconv = Dict{String,DataType}( #String to DataType conversion
    "Int32" => Int32, 
    "Int16" => Int16, 
    "Int8"  => Int8,
    "Float32" => Float32,
    "Float64" => Float64
)

"""
function get_varname(name::String,nc::NcFile)
Find name of a variable in an nc file using a list of alias values
Example: xname = get_name("x_node",keys(ncfile))
returns nothing if variable does not exist on the file
"""
function get_varname(name::String,nc::NcFile)
    ncvarnames=keys(nc)
    if !(name in keys(aliases))
        error("unknown variable")
    end
    tempstrs=intersect(aliases[name],ncvarnames)
    if length(tempstrs)==0
        return nothing
    else
        return first(tempstrs)
    end
end

"""
function get_coord_info(mapfiles::Vector{String})
Example: info=get_coord_info(["test_data/locxxz_map.nc"])
       ymin=info["ymin"]
Extract some info from mapfiles. The purpose is to collect the
info that is needed for the default config.  
"""
function get_coord_info(mapfiles::Vector{String})
    firstmap=NetCDF.open(first(mapfiles))
    # coordinate names
    xname_node=get_varname("x_node",firstmap)
    if xname_node===nothing
        error("Could not find x-coordinate variable in $(firstmap)")
    end
    yname_node=get_varname("y_node",firstmap)
    if yname_node===nothing
        error("Could not find y-coordinate variable in $(firstmap)")
    end
    xname_center=get_varname("x_center",firstmap)
    if xname_center===nothing
        error("Could not find x-coordinate variable in $(firstmap)")
    end
    yname_center=get_varname("y_center",firstmap)
    if yname_center===nothing
        error("Could not find y-coordinate variable in $(firstmap)")
    end
    zname_center=get_varname("z_center",firstmap)
    if zname_center===nothing
        zname_center=""
    end
    # determine bbox
    xs=firstmap[xname_node][:]
    xmax=maximum(xs)
    xmin=minimum(xs)
    ys=firstmap[yname_node][:]
    ymax=maximum(ys)
    ymin=minimum(ys)
    xs_center=firstmap[xname_center][:]
    ncells=[length(xs_center)]
    varnames=[]
    for varname in try_vars
        ncvar=get_varname(varname,firstmap)
        if !(ncvar===nothing)
            push!(varnames,varname)
        end
    end
    finalize(firstmap)
    # look at other map files
    for i=2:length(mapfiles)
        mapfile=NetCDF.open(mapfiles[i])
        xs=mapfile[xname_node][:]
        xmax=max(xmax,maximum(xs))
        xmin=min(xmin,minimum(xs))
        ys=mapfile[yname_node][:]
        ymax=max(ymax,maximum(ys))
        ymin=min(ymin,minimum(ys))
        xs_center=mapfile[xname_center][:]
        push!(ncells,length(xs_center))
        finalize(mapfile)
    end
    return Dict{String,Any}("xmin" => xmin, "ymin" => ymin, "xmax" => xmax, "ymax" => ymax,
                            "ncells" => ncells, "xname_node" => xname_node, "yname_node" => yname_node,
                            "xname_center" => xname_center, "yname_center" => yname_center, 
                            "zname_center" => zname_center, "varnames" => varnames)
end

"""
function default_config(mapfiles::Vector{String})
 Example: conf=default_config("myrun_map_interp.nc")
"""
function default_config(mapfiles::Vector{String})
    #firstmap=first(mapfiles)
    config=Dict{String,Any}()
    globals=Dict{String,Any}()
    globals["map_files"]=mapfiles
    zarrname=replace(first(mapfiles), ".nc" => ".zarr")
    zarrname=replace(zarrname, r"_[0-9]+" => s"")
    zarrname=replace(zarrname, r".*/" => s"")
    globals["zarr_file"]=zarrname
    globals["chunk_target_size"]=chunk_target_size
    info=get_coord_info(mapfiles)
    globals["extent"]=Dict{String,Any}("xmin" => info["xmin"], "xmax" =>info["xmax"],
                                       "ymin" => info["ymin"], "ymax" =>info["ymax"])
    globals["islatlon"]=false
    #propose dx and dy for regular grid 
    lx=info["xmax"]-info["xmin"]
    ly=info["ymax"]-info["ymin"]
    n=sum(info["ncells"])
    area_per_cell=lx*ly/n
    aspect_ratio=lx/ly
    dx=sqrt(area_per_cell) #estimate dy same
    nx=max(round(Int64,lx/dx)-1,1)
    ny=max(round(Int64,ly/dx)-1,1)
    globals["nx"]=nx
    globals["ny"]=ny
    config["global"]=globals
    varnames=info["varnames"]
    for varname in varnames
        varconfig=Dict{String,Any}(
            "scale_factor" => defaults[varname]["scale_factor"],
            "add_offset"   => defaults[varname]["add_offset"],
            "data_type"    => defaults[varname]["data_type"],
        )
        config[varname]=varconfig
   end
   return config
end

"""
function default_map_chunksize(nx,ny,ntarget)
 example: (ncx,ncy)=default_map_chunksize(nx,ny,ntarget)
 with 40 x 10 cells and a target size of 100, this will give 10x10 chunks
 """
function default_map_chunksize(nx,ny,ntarget)
    nx_chunks=max(1,round(Int64,sqrt(ntarget)))
    ny_chunks=max(1,round(Int64,sqrt(ntarget)))
    if nx_chunks>nx
        nx_chunks=nx
        ny_chunks=max(1,round(Int64,(ntarget/nx)))
        ny_chunks=min(ny_chunks)
    end
    if ny_chunks>ny
        ny_chunks=ny
        nx_chunks=max(1,round(Int64,(ntarget/ny)))
        nx_chunks=min(nx_chunks,nx)
    end
    return (nx_chunks,ny_chunks)
end

"""
function default_his_chunksize(Vector::varsize,ntarget=1000000)
    cz = default_his_chunksize((1000,1000),10000)
    Strategy is to make chunks over stations
"""
function default_his_chunksize(varsize,ntarget=1000000)
    blocksize=[varsize...] #make mutable
    rank=length(varsize)
    n=prod(varsize)
    blockpref=1
    if rank==3
       blockpref=2
    end
    if n>ntarget
       nblocks=div(n,ntarget)
       nblocks=max(1,nblocks)
       blocklen=div(varsize[blockpref],nblocks)
       blocklen=max(blocklen,1)
       blocksize[blockpref]=blocklen
    end
    return tuple(blocksize...)
 end
 
"""
function varlist(config::Dict{String,Any})
  vars=varlist(config)
"""
function varlist(config::Dict{String,Any})
    vars=Vector{String}()
    for varname in keys(config)
        if !startswith(lowercase(varname),r"global")
            push!(vars,varname)
        end
    end
    return vars
end

"""
function vardims(var::NcVar)
    example: mydims = vardims(nc.vars["waterlevel"])
"""
function vardims(var::NcVar)
    ndims=var.ndim
    return reverse([var.dim[i].name for i in 1:ndims])
end

"""
function copy_var(input::NcFile,output,varname,config,stoponmissing=true)
"""
function copy_var(input::NcFile,output,varname,config,stop_on_missing=true)
    println("copy variable name=$(varname)")
    ncname=get_varname(varname,input)
    if (ncname===nothing)
       if stop_on_missing
          error("could not find variable $(varname) in $(input.name).")
       else
          return nothing
       end
    end
    in_var=input.vars[ncname]
    in_atts=in_var.atts
    in_type=typeof(in_var).parameters[1]
    in_dummy=get(in_atts,"_FillValue",9999.0)
    in_size=size(in_var)
    in_rank=length(in_size)
    #output
    if !haskey(config,varname) #create empty one if absent
        config[varname]=Dict{String,Any}()
    end
    out_offset=get(config[varname],"add_offset",defaults[varname]["add_offset"])
    out_scale=get(config[varname],"scale_factor",defaults[varname]["scale_factor"])
    out_dummy=get(config[varname],"_FillValue",defaults[varname]["_FillValue"])
    out_type_str=get(config[varname],"data_type",defaults[varname]["data_type"])
    out_type=typeconv[out_type_str]
    out_chunk_target_size=get(config["global"],"chunk_target_size",chunk_target_size)
    out_chunk_size=default_his_chunksize(in_size,out_chunk_target_size)
    out_atts=copy(in_atts)
    out_atts["scale_factor"]=out_scale
    out_atts["add_offset"]=out_offset
    out_atts["_ARRAY_DIMENSIONS"]=vardims(in_var)
    ###DEBUG MVL 
    println("varname= $(varname)")
    println("in_dummy = $(in_dummy)")
    println("out_dummy= $(out_dummy)")
    #create output var
    out_var = zcreate(out_type, output, varname,in_size...,attrs=out_atts, chunks = out_chunk_size)
    println("in_size= $(in_size)")
    println("out_size= $(size(out_var))")
    #copy content
    if in_rank==1
        #in one go 
        in_temp=in_var[:]
        if out_type<:Int
            out_temp=round.(out_type,(in_temp.-out_offset)./out_scale)
            out_temp[in_temp.==in_dummy].=out_dummy
            out_var[:].=out_temp[:]
        else
            out_temp=(in_temp.-out_offset)./out_scale
            out_temp[in_temp.==in_dummy].=out_dummy
            out_var[:].=out_temp[:]    
        end
    elseif in_rank==2
        println("out_chunk_size = $(out_chunk_size)")
        if prod(in_size)==prod(out_chunk_size)
            #in one go 
            in_temp=in_var[:,:]
            out_temp=round.(out_type,(in_temp.-out_offset)./out_scale)
            out_temp[in_temp.==in_dummy].=out_dummy
            out_var[:,:].=out_temp[:,:]
        else #multiple blocks in time
            nblocks=max(div(prod(in_size),prod(out_chunk_size)),1)
            dimlen=in_size[2]
            blockstep=max(1,div(dimlen,nblocks))
            ifirst=1
            while ifirst<dimlen
                ilast=min(ifirst+blockstep-1,dimlen)
                in_temp=in_var[:,ifirst:ilast]
                out_temp=round.(out_type,(in_temp.-out_offset)./out_scale)
                out_temp[in_temp.==in_dummy].=out_dummy
                out_var[:,ifirst:ilast]=out_temp[:,:]
            end
        end
    elseif in_rank==3
        if prod(in_size)==prod(out_chunk_size)
            #in one go 
            in_temp=in_var[:,:,:]
            out_temp=round.(out_type,(in_temp.-out_offset)./out_scale)
            out_temp[in_temp.==in_dummy].=out_dummy
            out_var[:,:,:].=out_temp[:,:,:]
        else #loop over multiple blocks in time, even though output has different chunks
            nblocks=max(div(prod(in_size),prod(out_chunk_size)),1)
            dimlen=in_size[3]
            blockstep=max(1,div(dimlen,nblocks))
            ifirst=1
            while ifirst<dimlen
                ilast=min(ifirst+blockstep-1,dimlen)
                in_temp=in_var[:,:,ifirst:ilast]
                out_temp=round.(out_type,(in_temp.-out_offset)./out_scale)
                out_temp[in_temp.==in_dummy].=out_dummy
                out_var[:,:,ifirst:ilast]=out_temp[:,:]
            end
        end
    else
        error("Wrong rank ($(in_rank)) for variable $(varname)")
    end
end

"""
function scale_values(in_values,in_dummy,out_type,out_offset,out_scale,out_dummy)
Convert array from float types to integer type
Values in teh input that are NaN or equal to in_dummy are set to out_dummy 
"""
# function scale_values(in_values,in_dummy,out_type,out_offset,out_scale,out_dummy)
# # This implementation claims too much memory
#    in_dummies=in_values.==in_dummy
#    in_nans=isnan.(in_values)
#    in_values[in_dummies].=out_offset
#    in_values[in_nans].=out_offset
#    out_max=typemax(out_type)
#    out_min=typemin(out_type)
#    temp=(in_values.-out_offset)./out_scale
#    temp=min.(temp,out_max)
#    temp=max.(temp,out_min)
#    out_temp=round.(out_type,temp)
#    out_temp[in_dummies].=out_dummy
#    out_temp[in_nans].=out_dummy
#    return out_temp
# end

function scale_values(in_values,in_dummy,out_type,out_offset,out_scale,out_dummy)
    # optimized version
    out_max=typemax(out_type)
    out_min=typemin(out_type)
    out_values=Array{out_type}(undef,size(in_values))
    for i in eachindex(in_values)
        in_value = in_values[i]
        if isnan(in_value) || (in_value==in_dummy)
            out_value = out_dummy
        else
            temp_value = (in_value - out_offset)/out_scale
            temp_value = min(temp_value,out_max)
            temp_value = max(temp_value,out_min)
            out_value = round(out_type,temp_value)
        end
        out_values[i]=out_value
    end
    return out_values
 end


function interp_var(inputs::Vector{NcFile},interp::Interpolator,output::ZGroup,varname::String,xpoints,ypoints,config,dumval=NaN)
    println("interpolating variable name=$(varname)")
    globals=config["global"]
    nx=length(xpoints)
    ny=length(ypoints)
    if !haskey(config,varname) #create empty one if absent
        config[varname]=Dict{String,Any}()
    end
    #output time and scaling
    out_offset=get(config[varname],"add_offset",defaults[varname]["add_offset"])
    out_scale=get(config[varname],"scale_factor",defaults[varname]["scale_factor"])
    out_dummy=get(config[varname],"_FillValue",defaults[varname]["_FillValue"])
    out_type_str=get(config[varname],"data_type",defaults[varname]["data_type"])
    out_type=typeconv[out_type_str]
    println("  write output as $(out_type_str) scaling $(out_scale) offset $(out_offset)")
    # compute size of chunks
    chunk_target_size=globals["chunk_target_size"]
    (x_chunksize,y_chunksize)=default_map_chunksize(nx,ny,chunk_target_size)
    println("(nx,ny)=$((nx,ny))  chunks=$((x_chunksize,y_chunksize))")
    # name of variable on netcdf file
    firstmap=first(inputs)
    nt=length(firstmap["time"])
    ncname=get_varname(varname,firstmap)
    if ncname===nothing
        error("could not find variable $(varname) in $(firstmap.name).")
    end
    in_dims=vardims(firstmap.vars[ncname])
    hastime="time" in in_dims
    is2d=(length(in_dims)==2)
    if !hastime
        is2d=(length(in_dims)==1)
    end
    in_size=size(firstmap.vars[ncname])
    # create variable and copy attributes 
    varatts=firstmap[ncname].atts
    varatts["scale_factor"]=out_scale
    varatts["add_offset"]=out_offset
    in_dummy=get(varatts,"_FillValue",9999.0)
    varatts["_FillValue"]=out_dummy
    #println("  attributes = $(varatts)")
    
    if is2d
        if hastime
            varatts["_ARRAY_DIMENSIONS"]=["time","y","x"]
            varatts["coordinates"]="time y_center x_center"
            var = zcreate(out_type, output, varname,(nx,ny,nt)...,attrs=varatts,chunks=(x_chunksize,y_chunksize,1))
            print("times $(nt):")
            for it=1:nt
                print("|")
                in_temp_uninterpolated=load_nc_map_slice(inputs,ncname,it)
                in_temp=interpolate(interp,xpoints,ypoints,in_temp_uninterpolated,dumval)
                out_temp=scale_values(in_temp,in_dummy,out_type,out_offset,out_scale,out_dummy)
                var[:,:,it]=out_temp[:,:]
            end
        else
            varatts["_ARRAY_DIMENSIONS"]=["y","x"]
            varatts["coordinates"]="y_center x_center"
            var = zcreate(out_type, output, varname,(nx,ny)...,attrs=varatts,chunks=(x_chunksize,y_chunksize))
            print("time independent")
            it=0
            in_temp_uninterpolated=load_nc_map_slice(inputs,ncname,it)
            in_temp=interpolate(interp,xpoints,ypoints,in_temp_uninterpolated,dumval)
            out_temp=scale_values(in_temp,in_dummy,out_type,out_offset,out_scale,out_dummy)
            var=out_temp
        end
    else #is 3D
        if hastime #x y z t
            nz=in_size[1]
            varatts["_ARRAY_DIMENSIONS"]=["time","z","y","x"]
            varatts["coordinates"]="time z_center_3d y_center x_center"
            at_center=true
            if varname in znames_iface
                varatts["_ARRAY_DIMENSIONS"]=["time","z_iface","y","x"]
                varatts["coordinates"]="time z_iface_3d y_center x_center"
                at_center=false
            end
            if (nx*ny*nz)<chunk_target_size
                var = zcreate(out_type, output, varname,(nx,ny,nz,nt)...,attrs=varatts,chunks=(x_chunksize,y_chunksize,nz,1))
            else
                var = zcreate(out_type, output, varname,(nx,ny,nz,nt)...,attrs=varatts,chunks=(x_chunksize,y_chunksize,1,1))
            end
            print("times $(nt):")
            for it=1:nt
                print("|")
                for ilayer in 1:nz
                    print(".")
                    #in_temp_uninterpolated=[]
                    # check if variable is given on edges
                    if varatts["location"]=="edge" #not face
                        @time in_temp_uninterpolated=load_nc_map_slice_at_faces(inputs,ncname,it,ilayer)
                    else
                        @time in_temp_uninterpolated=load_nc_map_slice(inputs,ncname,it,ilayer)
                    end
                    @time in_temp=interpolate(interp,xpoints,ypoints,in_temp_uninterpolated,dumval)
                    @time out_temp=scale_values(in_temp,in_dummy,out_type,out_offset,out_scale,out_dummy)
                    var[:,:,ilayer,it]=out_temp[:,:]
                end
            end
        else # x y z
            nz=in_size[1]
            varatts["_ARRAY_DIMENSIONS"]=["z","y","x"]
            varatts["coordinates"]="z_center y_center x_center"
            if varname in znames_iface
                varatts["_ARRAY_DIMENSIONS"]=["z_iface","y","x"]
                varatts["coordinates"]="z_iface y_center x_center"
            end
            if (nx*ny*nz)<chunk_target_size
                var = zcreate(out_type, output, varname,(nx,ny,nz)...,attrs=varatts,chunks=(x_chunksize,y_chunksize,nz,1))
            else
                var = zcreate(out_type, output, varname,(nx,ny,nz)...,attrs=varatts,chunks=(x_chunksize,y_chunksize,1))
            end
            print("time independent")
            for ilayer in 1:nz
                print(".")
                it=0
                in_temp_uninterpolated=load_nc_map_slice(inputs,ncname,it,ilayer)
                in_temp=interpolate(interp,xpoints,ypoints,in_temp_uninterpolated,dumval)
                out_temp=scale_values(in_temp,in_dummy,out_type,out_offset,out_scale,out_dummy)
                var[:,:,ilayer]=out_temp[:,:]
            end
        end
    end
    println("")
end


"""
function copy_strings(input::NcFile,output,varname,config)
"""
function copy_strings(input::NcFile,output,varname,config)
    #input
    in_var=input.vars[varname]
    in_atts=in_var.atts
    in_type=typeof(in_var).parameters[1]
    if in_type!=NetCDF.ASCIIChar
        error("Can only copy character arrays to Strings, varname = $(varname)")
    end
    in_size=size(in_var)
    in_rank=length(in_size)
    if in_rank!=2
        error("Only 2d character arrays are supported for now, varname = $(varname)")
    end
    in_stations=in_size[2]
    in_strlen=in_size[1]
    println("in_stations= $(in_stations)")
    #output
    out_type=Zarr.MaxLengthString{in_strlen,UInt8}
    out_atts=copy(in_atts)
    out_atts["_ARRAY_DIMENSIONS"]=[vardims(in_var)[1]]
    #create output ar
    out_var = zcreate(out_type, output, varname,in_stations,attrs=out_atts)
    #copy content
    #in one go
    in_temp=in_var[:,:] 
    out_temp=replace.(nc_char2string(in_temp),r"\s+" => "")
    println("out_temp= $(out_temp)")
    out_var[:].=out_temp[:]    
end

"""
function has_z(input::NcFile,vars::Vector)
Report if any of the variables has a z coordinate.
We now take a shortcut and check if the variable is a 3d array.
"""
function has_z(input::NcFile,vars::Vector)
    result=false
    for varname in vars
        r=length(size(input.vars[varname]))
        if r==3
            result=true
        end
    end
    return result
end

#
# main function
#
function main(args)
    # read user input
    n = length(args)
    if n < 1
        println("Usage: julia dflow_map_interp_to_zarr.jl config.toml")
        println("This command will generate a config file:")
        println("julia dflow_map_interp_to_zarr.jl myrun_map_00*.nc")
    elseif endswith(lowercase(first(args)),r".toml")
        configfile=first(args)
        if !isfile(configfile)
            throw(ArgumentError("File not found: $(configfile)\n"))
        end
        # load configfile
        config=TOML.parsefile(configfile)
        if !haskey(config,"global")
            error("Missing [global] group in config file")
        end
        globals=config["global"]
        if !haskey(globals,"map_files")
            error("Missing keyword map_files in group [global] in config file")
        end 
        map_files=globals["map_files"]
        if !isfile(first(map_files))
            error("Cannot find file: $(first(map_files))")
        end
        outname=""
        if haskey(globals,"zarr_file")
            outname=globals["zarr_file"]
        else
            error("missing key zarr_file in configuration)")
        end
        if ispath(outname)
            error("Output name $(outname) exists. Will not overwrite.")
        end
        # initialize netcdf maps files
        map=[NetCDF.open(filename) for filename in map_files]
        firstmap=first(map)
        # Create zarr file and start copying metadata
        globalattrs = firstmap.gatts
        output = zgroup(outname,attrs=globalattrs)
        # output coordinates
        islatlon=globals["islatlon"]
        nx=max(1,globals["nx"]-1)
        ny=max(1,globals["ny"]-1)
        extent=globals["extent"]
        xmin=extent["xmin"]
        xmax=extent["xmax"]
        ymin=extent["ymin"]
        ymax=extent["ymax"]
        dx=(xmax-xmin)/nx
        dy=(ymax-ymin)/ny
        xmin=xmin+0.5*dx
        xmax=xmax-0.5*dx
        ymin=ymin+0.5*dy
        ymax=ymax-0.5*dy
        xpoints=collect(range(xmin,stop=xmax,length=nx))
        ypoints=collect(range(ymin,stop=ymax,length=ny))
        println("x: $(xmin) : $(dx) $(xmax) : $(nx)")
        println("y: $(ymin) : $(dy) $(ymax) : $(ny)")
        # write coordinates
        xname_node=get_varname("x_node",firstmap)
        xatts=firstmap[xname_node].atts 
        xatts["_ARRAY_DIMENSIONS"]=["x"]
        xvar = zcreate(Float64, output, "x_center",length(xpoints),attrs=xatts)
        xvar[:]=xpoints[:]
        yname_node=get_varname("y_node",firstmap)
        yatts=firstmap[yname_node].atts 
        yatts["_ARRAY_DIMENSIONS"]=["y"]
        yvar = zcreate(Float64, output, "y_center",length(ypoints),attrs=yatts)
        yvar[:]=ypoints[:]
        # set up interpolation
        interp=load_dflow_grid(map,50,islatlon)
        # interpolate per variable
        vars=varlist(config)
        for varname in vars
            println("Interpolating for variable $(varname)")
            @time interp_var(map,interp,output,varname,xpoints,ypoints,config) #TODO Far too much memory is allocated here!
        end
        #copy dimensions and coordinates
        # make additional fullgrid z coordinates options in try_vars above
        # nc_zc=get_varname("z_center_3d",firstmap) #try 3d z coords
        # if !(nc_zc==nothing)
        #     interp_var(map,interp,output,"z_center_3d",xpoints,ypoints,config)
        # end
        # nc_zw=get_varname("z_iface_3d",firstmap) #try 3d z coords
        # if !(nc_zw==nothing)
        #     interp_var(map,interp,output,"z_iface_3d",xpoints,ypoints,config)
        # end
        copy_var(firstmap,output,"z_center",config,false) #try 1d z coords
        copy_var(firstmap,output,"time",config)
        # create consolidate_metadata for faster internet access
        Zarr.consolidate_metadata(output)
        
    else # expect list of mapfiles and generate default config
        first_mapfile=first(args)
        if !isfile(first_mapfile)
            throw(ArgumentError("File not found: $(first_mapfile)\n"))
        end
        if !endswith(lowercase(first(args)),r"_map.nc")
            throw(ArgumentError("Expecting a Deflt3D-FM map netcdf file: $(first_mapfile)\n"))
        end
        configfile="config_maps_interp.toml"
        if isfile(configfile)
            throw(ArgumentError("Existing config file. Will not overwrite: $(configfile)\n"))
        end
        config=default_config(args)
        configfile="config_maps_interp.toml"
        open(configfile, "w") do io
            TOML.print(io, config)
        end
        if(debuglevel>0)
            TOML.print(config)
        end
        @info "Configuration has been written to $(configfile)"
    end
end


#
# Call main if used as a script, but not if loaded as a module 
#

# some defaults for manual tesing
mapfiles=["test_data/estuary_0000_map.nc", "test_data/estuary_0000_map.nc"]
#mapfiles=["test_data/locxxz_map.nc"]
configfile=["config_map_interp.toml"] #TODO these names are not used

# do nothing when called as module
if abspath(PROGRAM_FILE) == @__FILE__
    println("ARGS = $(ARGS)")
    @time main(ARGS)
end

nothing
