# This script downloads the relevant data for this case. This data is not inside the repository 
# because of its size.

using HTTP
#170Mb wave data from SWAN model
# curl -o swan_20171007_28.nc  https://nx7384.your-storageshare.de/s/WPcFEfgtqA5bnsS/download
HTTP.download("https://nx7384.your-storageshare.de/s/WPcFEfgtqA5bnsS/download","swan_20171007_28.nc")
