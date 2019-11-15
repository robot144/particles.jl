# This script downloads the relevant data for this case. This data is not inside the repository 
# because of its size.

using HTTP

# Winds from ERA-5  for this area and timespan
# 165Mb
HTTP.download("https://nx7384.your-storageshare.de/s/HosJSLErsgLd4TT/download","cera5_wind_201403_05.nc") 

# Tidal currents from GTSM3 with 5km coastal resolution
# only domain 0001 of size 15Gb
HTTP.download("https://nx7384.your-storageshare.de/s/xFRpoFJMXzNeT96/download","gtsm_fine_0001_map.nc") 
