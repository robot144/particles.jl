#! /bin/sh
# This script downloads the relevant data for this case. This data is not inside the repository 
# because of its size.

# Winds from ERA-5  for this area and timespan
# 165Mb
curl -o era5_wind_201403_05.nc https://nx7384.your-storageshare.de/s/HosJSLErsgLd4TT/download

# Tidal currents from GTSM3 with 5km coastal resolution
# only domain 0001 of size 15Gb
curl -o gtsm_fine_0001_map.nc https://nx7384.your-storageshare.de/s/xFRpoFJMXzNeT96/download

# Ocean currents and sea-ice coverage and velocities from CMEMS global reanalysis
# size 138Mb
curl -o cmems_201403_05.nc  https://nx7384.your-storageshare.de/s/29NxGwBeS46xe5T
