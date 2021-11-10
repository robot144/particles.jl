#! /bin/bash
# This script downloads the relevant data for this case. This data is not inside the repository 
# because of its size.

# CMEMS flow field
# 223 Mb per file 
wget -O u0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc https://nx7384.your-storageshare.de/s/BG3reLnBJktAMwz/download &
wget -O v0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc https://nx7384.your-storageshare.de/s/KgKsdmWk5FqbZPP/download

# winds from NOAA GFS
# 48Mb
wget -O gfs_winds.nc https://nx7384.your-storageshare.de/s/K34zDXKK2GpPLBM/download
