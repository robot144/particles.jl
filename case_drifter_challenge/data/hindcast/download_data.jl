# This script downloads the relevant data for this case. This data is not inside the repository 
# because of its size.

using HTTP

# CMEMS flow field
# 223 Mb per file 
HTTP.download("https://nx7384.your-storageshare.de/s/BG3reLnBJktAMwz/download","u0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/KgKsdmWk5FqbZPP/download","v0_2021-05-29_00-00-00_2021-06-22_00-00-00.nc")

# winds from NOAA GFS
# 48Mb one day
#HTTP.download("https://nx7384.your-storageshare.de/s/K34zDXKK2GpPLBM/download","gfs_winds.nc")
# 763Mb 31 days
HTTP.download("https://nx7384.your-storageshare.de/s/jMQtERxRakj59ii/download,"gfs_winds.nc")
