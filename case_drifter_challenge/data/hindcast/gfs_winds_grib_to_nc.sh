#! /bin/bash
# Convert multiple grib files to a single netcdf, keeping only the variables 10u aud 10v
#
# makes use of ecmwf's eccodes. On ubuntu install with "sudo apt install libeccodes-tools"
# and cdo. On ubuntu install with "sudo apt install cdo"
#

# first create temp folder or stop
if [ ! -d temp ]; then
	echo "creating folder temp"
	mkdir temp
else
	echo "ERROR: folder temp exists. Remove if you want to continue."
	exit 1
fi

# select 10u and 10v
for file in *.grb2 ;do
	ufile="${file/.grb2/_10u.grb2}"
	grib_copy -w shortName=10u "$file" "temp/$ufile"
	vfile="${file/.grb2/_10v.grb2}"
	grib_copy -w shortName=10v "$file" "temp/$vfile"
done

#convert to netcdf
#grib_to_netcdf -o gfs_winds_10u.nc temp/*_10u.grb2 #PROBLEMATIC
#grib_to_netcdf -o gfs_winds_10v.nc temp/*_10v.grb2
cdo -f nc4 mergetime temp/*_10u.grb2 gfs_winds_10u.nc
cdo -f nc4 mergetime temp/*_10v.grb2 gfs_winds_10v.nc

cdo merge gfs_winds_10u.nc gfs_winds_10v.nc gfs_winds.nc
