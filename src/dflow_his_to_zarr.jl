#
# Convert delft3d time-series to zarr format
# This add selection, compression and chuning into separate files per chunk
#
using TOML
using NetCDF
using Zarr
using Dates

debuglevel=1

#
# defaults
#
try_vars = ["waterlevel","salinity","temperature","x_velocity","y_velocity","vicww"] #variables to look for in hisfile
defaults = Dict(
    "waterlevel" => Dict(
       "scale_factor" => 0.001, 
       "add_offset" => 0.0,
       "data_type" => "Int16",
       "_FillValue" => 9999 ),
    "salinity" => Dict(
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "temperature" => Dict(
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "x_velocity" => Dict(
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "y_velocity" => Dict(
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "vicww" => Dict(
        "scale_factor" => 1e-6, 
        "add_offset" => 0.0,
        "data_type" => "Int16",
        "_FillValue" => 9999 ),
    "time" => Dict(
        "scale_factor" => 1.0, 
        "add_offset" => 0.0,
        "data_type" => "Float64",
        "_FillValue" => -9999.0),
    "station_x_coordinate" => Dict(
        "scale_factor" => 1.0, 
        "add_offset" => 0.0,
        "data_type" => "Float64",
        "_FillValue" => -9999.0),
    "station_y_coordinate" => Dict(
        "scale_factor" => 1.0, 
        "add_offset" => 0.0,
        "data_type" => "Float64",
        "_FillValue" => -9999.0),
    "zcoordinate_c" => Dict(
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int32",
        "_FillValue" => -9999.0),
    "zcoordinate_w" => Dict(
        "scale_factor" => 0.01, 
        "add_offset" => 0.0,
        "data_type" => "Int32",
        "_FillValue" => -9999.0)
    )
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
function default_config(hisfile::String)
 Example: conf=default_config("myrun_his.nc")
"""
function default_config(hisfile::String)
    config=Dict{String,Any}()
    globals=Dict{String,Any}()
    globals["history_file"]=hisfile
    zarrname=replace(hisfile, ".nc" => ".zarr")
    zarrname=replace(zarrname, r".*/" => s"")
    globals["zarr_file"]=zarrname
    globals["chunk_target_size"]=chunk_target_size
    config["global"]=globals
    nc=NetCDF.open(hisfile)
    ncvars=keys(nc)
    for varname in try_vars
        if varname in ncvars
            varconfig=Dict{String,Any}(
                "scale_factor" => defaults[varname]["scale_factor"],
                "add_offset"   => defaults[varname]["add_offset"],
                "data_type"    => defaults[varname]["data_type"],
            )
            config[varname]=varconfig
        end
   end
   finalize(nc)
   return config
end

"""
function range_in_chunks(total_range::OrdinalRange, chunksize::Int)
Divide a range line 1:2:10 into chunks, eg [1:2:5,7:2:9]
This is often usefull to make a loop where the computations are performed per chunk.
"""
function range_in_chunks(total_range::OrdinalRange, chunksize::Int)
    chunks=Vector{StepRange}()
    n=length(total_range)
    n_chunks=max(ceil(Int64,n/chunksize),1)
    println("#chunks $(n_chunks) chunk size $(chunksize)")
    for i=1:n_chunks
       istart=total_range[(i-1)*chunksize+1]
       istep=((total_range isa UnitRange) ? 1 : total_range.step)
       if i<n_chunks
          istop=total_range[i*chunksize]
       else
          istop=total_range[end]
       end
       chunk=istart:istep:istop
       push!(chunks,chunk)
    end
    return chunks
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
            in_temp[in_temp.==in_dummy].=out_offset
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
            in_temp[in_temp.==in_dummy].=out_offset
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
                in_temp[in_temp.==in_dummy].=out_offset
                out_temp=round.(out_type,(in_temp.-out_offset)./out_scale)
                out_temp[in_temp.==in_dummy].=out_dummy
                out_var[:,ifirst:ilast]=out_temp[:,:]
                ifirst=ilast
            end
        end
    elseif in_rank==3
        if prod(in_size)==prod(out_chunk_size)
            #in one go 
            in_temp=in_var[:,:,:]
            in_temp[in_temp.==in_dummy].=out_offset
            out_temp=round.(out_type,(in_temp.-out_offset)./out_scale)
            out_temp[in_temp.==in_dummy].=out_dummy
            out_var[:,:,:].=out_temp[:,:,:]
        else #loop over multiple blocks in time, even though output has different chunks
            #buffer_target_size=10^8
            buffer_target_size=10^9
            n_per_time=prod(in_size[1:2])
            blockstep=max(1,ceil(Int64,buffer_target_size/n_per_time))
            dimlen=in_size[3]
            for chunk in range_in_chunks(1:dimlen,blockstep)
                ifirst=chunk.start
                ilast=chunk.stop
                thislen=ilast-ifirst+1
                println("   copying times $(ifirst) - $(ilast) of $(dimlen)")
                print("<R")
                in_temp=zeros(in_size[1],in_size[2],thislen)
                # in_temp=in_var[:,:,ifirst:ilast] # all in one 
                for i in chunk
                   in_temp[:,:,(i-ifirst+1)]=in_var[:,:,i]
                   if rem(i,100)==0
                      print(".")
                   end
                end
                print(">")
                bool_dummy=in_temp.==in_dummy
                bool_nan=isnan.(in_temp)
                in_temp[bool_dummy].=out_offset
                in_temp[bool_nan].=out_offset
                out_temp=round.(out_type,(in_temp.-out_offset)./out_scale)
                out_temp[bool_dummy].=out_dummy
                out_temp[bool_nan].=out_dummy
                print("<W")
                out_var[:,:,ifirst:ilast]=out_temp[:,:,:]
                println(">")
                ifirst=ilast
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
        println("Usage: julia dflow_his_to_zarr.jl config.toml")
        println("This command will generate a config file:")
        println("julia dflow_his_to_zarr.jl myrun_his.nc")
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
        if !haskey(globals,"history_file")
            error("Missing keyword history_file in group [global] in config file")
        end 
        hisfile=globals["history_file"]
        if !isfile(hisfile)
            error("Cannot find file: $(hisfile)")
        end
        his=NetCDF.open(hisfile)
        vars=varlist(config)
        if !endswith(lowercase(hisfile),r"_his.nc")
            error("History file $(hisfile) is not like _his.nc")
        end
        outname=replace(hisfile, ".nc" => ".zarr")
        if haskey(globals,"zarr_file")
            outname=globals["zarr_file"]
        end
        if ispath(outname)
            error("Output name $(outname) exists. Will not overwrite.")
        end
        globalattrs = his.gatts
        out = zgroup(outname,attrs=globalattrs)
        for varname in vars
            println("copying variable $(varname)")
            copy_var(his,out,varname,config)
        end
        #copy dimensions and coordinates
        copy_var(his,out,"time",config)
        copy_var(his,out,"station_x_coordinate",config)
        copy_var(his,out,"station_y_coordinate",config)
        copy_strings(his,out,"station_name",config)
        if has_z(his,vars)
            copy_var(his,out,"zcoordinate_c",config)
            copy_var(his,out,"zcoordinate_w",config)
        end
        Zarr.consolidate_metadata(out)
        ## TODO 
    else # expect history file and generate default config
        hisfile=first(args)
        if !isfile(hisfile)
            throw(ArgumentError("File not found: $(hisfile)\n"))
        end
        if !endswith(lowercase(first(args)),r"_his.nc")
            throw(ArgumentError("Expecting a Deflt3D-FM history netcdf file: $(hisfile)\n"))
        end
        configfile="config_his.toml"
        if isfile(configfile)
            throw(ArgumentError("Existing config file. Will not overwrite: $(configfile)\n"))
        end
        config=default_config(hisfile)
        configfile="config_his.toml"
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
hisfile=["test_data/locxxz_his.nc"]
configfile=["config_his.toml"]
#main(ARGS)

nothing
