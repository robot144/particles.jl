# IntZarrct with dflow netcdf output map files
#
using Dates
using ZipFile

#
# support functions
#

# function unzip(file,exdir="")
# extract all files from zip file to exdir
function unzip(file,exdir="")
    fileFullPath = isabspath(file) ?  file : joinpath(pwd(),file)
    basePath = dirname(fileFullPath)
    outPath = (exdir == "" ? basePath : (isabspath(exdir) ? exdir : joinpath(pwd(),exdir)))
    isdir(outPath) ? "" : mkdir(outPath)
    zarchive = ZipFile.Reader(fileFullPath)
    for f in zarchive.files
        fullFilePath = joinpath(outPath,f.name)
        if (endswith(f.name,"/") || endswith(f.name,"\\"))
            mkdir(fullFilePath)
        else
            write(fullFilePath, read(f))
        end
    end
    close(zarchive)
end


#
# test
#

function test1()
    #first unzip test data
    if isfile("../test_data/locxxz_fullgrid_map.zip")
        if !isdir("./locxxz_fullgrid_map.zarr")
            unzip("../test_data/locxxz_fullgrid_map.zip",".")
        end
    else
        error("test data not found in ../test_data/locxxz_fullgrid_map.zip")
    end

    Zarr_data = ZarrData(".","locxxz_fullgrid_map.zarr")

    p1=load_map_slice(Zarr_data,"waterlevel",2)
    @test size(p1)==(344,1)
    @test p1[1,1]≈-0.049

    t0=get_reftime(Zarr_data)
    @test t0==DateTime(2001,1,1)

    times=get_times(Zarr_data,t0)
    @test times[1]≈0.0
    @test times[2]≈20.1
    @test times[3]≈40.2 #rounding error in time array with long interval in seconds

    t2=as_DateTime(Zarr_data,t0,times[2])
    @test t2==DateTime("2001-01-01T00:00:20.1")

    h=initialize_interpolation(Zarr_data,"waterlevel",t0)
    h1=h(140.0,0.5,0.0,20.0) #x,y,z,t
    @test h1≈-0.0248756218
    h2=h(140.0,0.5,0.0,120.0)
    @test h2≈ 0.0079666666

    # TODO: 3d interpolation for z too
    u=initialize_interpolation(Zarr_data,"x_velocity",t0)
    u1=u(140.0,0.5,-2.0,20.0) #x,y,z,t
    @test u1≈0.024742951907131133
    u2=u(140.0,0.5,-2.0,120.0) #x,y,z,t
    @test u2≈-0.52095
end

test1()
