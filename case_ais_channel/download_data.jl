# This script downloads the relevant data for this case. This data is not inside the repository 
# because of its size.

using HTTP
# GTSM-FM curents 2D
# 1.6 Gb 
# curl -o gtsm_fine_0002_map.nc https://nx7384.your-storageshare.de/s/Qtw2jrc9PP243a8/download 
HTTP.download("https://nx7384.your-storageshare.de/s/Qtw2jrc9PP243a8/download","gtsm_fine_0002_map.nc")
