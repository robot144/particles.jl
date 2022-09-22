#
# Convert delft3d time-series to zarr format
# This add selection, compression and chuning into separate files per chunk
#
using TOML
using NetCDF
using Zarr
using Dates
include(joinpath(@__DIR__,"unstructured_grid.jl")) #bit ad-hoc
include(joinpath(@__DIR__,"dflow.jl")) #bit ad-hoc

debuglevel=1

#
# defaults
#
try_vars = ["waterlevel","x_velocity","y_velocity","salinity"] #variables to look for in hisfile
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
    "salinity" => Dict(
        "name" => "salinity",
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "time" => Dict(
        "scale_factor" => 1.0, 
        "add_offset" => 0.0,
        "data_type" => "Float64",
        "_FillValue" => -9999.0)
        )

# variables appear under different names in the delft3d-fm output files. Here we list the options
aliases=Dict{String,Vector{String}}(
    "waterlevel" => ["s0", "mesh2d_s0"],
    "x_velocity" => ["ucx", "mesh2d_ucx"],
    "y_velocity" => ["ucy", "mesh2d_ucy"],
    "salinity"   => ["sal","mesh2d_sal"],
    "x_center"   => ["FlowElem_xcc","mesh2d_face_x"],
    "y_center"   => ["FlowElem_ycc","mesh2d_face_y"],
    "z_center"   => ["mesh2d_layer_z"],
    "x_node"     => ["mesh2d_node_x","NetNode_x"],
    "y_node"     => ["mesh2d_node_y","NetNode_y"]
)
#znames_iface=["mesh2d_layer_z"]

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
function get_varname(name::String)
Find name of a variable in an nc file using a list of alias values
Example: xname = get_name("x_node",keys(ncfile))
returns nothing if variable does not exist on the file
"""
function get_varname(name::String,ncvarnames)
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
    varnames=keys(firstmap)
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
function default_config(hisfile::String)
 Example: conf=default_config("myrun_his.nc")
"""
function default_config(mapfiles::Vector{String})
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
    ncvarnames=info["varnames"]
    for varname in try_vars
        ncvarname=get_varname(varname,ncvarnames)
        if !(ncvarname===nothing)
           varconfig=Dict{String,Any}(
                "scale_factor" => defaults[varname]["scale_factor"],
                "add_offset"   => defaults[varname]["add_offset"],
                "data_type"    => defaults[varname]["data_type"],
            )
            config[varname]=varconfig
        end
   end
   return config
end

function default_map_nchunks(nx,ny,ntarget)
    nx_chunks=max(1,round(Int64,nx/sqrt(ntarget))) #product of chuns in both directions
    ny_chunks=max(1,round(Int64,ny/sqrt(ntarget)))
    if nx==1
        ny_chunks=max(1,round(Int64,ny/ntarget))
    elseif ny==1
        nx_chunks=max(1,round(Int64,nx/ntarget))
    end
    return (nx_chunks,ny_chunks)
end

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
function copy_var(input::NcFile,output,varname,config)
"""
function copy_var(input::NcFile,output,varname,config)
    #input
    in_var=input.vars[varname]
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
    #create output ar
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
        map=[NetCDF.open(filename) for filename in mapfiles]
        firstmap=first(map)
        # Create zarr file and start copying metadata
        globalattrs = firstmap.gatts
        out = zgroup(outname,attrs=globalattrs)
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
        # compute size of chunks
        chunk_target_size=globals["chunk_target_size"]
        (x_chunksize,y_chunksize)=default_map_chunksize(nx,ny,chunk_target_size)
        # write coordinates
        xname_node=get_varname("x_node",firstmap)
        xatts=firstmap[xname_node].atts 
        println("xatts = $(xatts)")
        xvar = zcreate(Float64, out, "x_center",length(xpoints),attrs=xatts)
        xvar[:]=xpoints[:]
        yname_node=get_varname("y_node",firstmap)
        yatts=firstmap[yname_node].atts 
        println("yatts = $(yatts)")
        yvar = zcreate(Float64, out, "y_center",length(ypoints),attrs=yatts)
        yvar[:]=ypoints[:]
        # set up interpolation
        interp=load_dflow_grid(map,50,islatlon)

        # interpolate per variable
        nt=length(firstmap["time"])
        vars=varlist(config)
        for varname in vars
            ncname=get_varname(varname,firstmap)
            varatts=firstmap[ncname].atts
            println("(nx,ny)=$((nx,ny))  chunks=$((x_chunksize,y_chunksize))")
            var = zcreate(Float64, out, varname,(nx,ny,nt)...,attrs=varatts,chunks=(x_chunksize,y_chunksize,1))
            for it=1:nt
                temp=load_nc_map_slice(map,ncname,it)
                temp_interp=interpolate(interp,xpoints,ypoints,temp)
                var[:,:,it]=temp_interp[:,:]
            end
        end
        #copy dimensions and coordinates
        copy_var(firstmap,out,"time",config)
        #TODO
        #if has_z(his,vars)
        #    copy_var(his,out,"zcoordinate_c",config)
        #    copy_var(his,out,"zcoordinate_w",config)
        #end
        Zarr.consolidate_metadata(out)
        
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
# main 
#
mapfiles=["test_data/estuary_0000_map.nc", "test_data/estuary_0000_map.nc"]
#mapfiles=["test_data/locxxz_map.nc"]
configfile=["config_map_interp.toml"]
#main(ARGS)

nothing
