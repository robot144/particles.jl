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

    p1=load_map_slice(Zarr_data,"waterlevel",1)
    @test size(p1)==(344,1)
    #@test p1[1,1]≈102574.4477038286

    # t0=get_reftime(Zarr_data)
    # @test t0==DateTime(2014,3,1)

    # times=get_times(Zarr_data,t0)
    # @test times[1]≈0.0
    # @test times[2]≈3600.0

    # t2=as_DateTime(Zarr_data,t0,times[2])
    # @test t2==DateTime("2014-03-01T01:00:00")

    # # u,v interpolation functions
    # t0=get_reftime(Zarr_data)
    # @test t0==DateTime(2014,3,1)

    # p=initialize_interpolation(Zarr_data,"msl",t0)
    # p1=p(20.0,75.0,0.0,0.0)
    # @test p1≈101877.70372941876
    # p2=p(20.0,75.0,0.0,3600.0)
    # @test p2≈101926.02316008728
end


test1()
