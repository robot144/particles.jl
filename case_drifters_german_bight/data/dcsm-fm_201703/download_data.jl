# This script downloads the relevant data for this case. This data is not inside the repository 
# because of its size.

using HTTP

# DCSM-FM-3D for March-April 2017
# 300Gb about 30Gb per domain
HTTP.download("https://nx7384.your-storageshare.de/s/2gFqaqsAwswic54/download","DCSM-FM_05nm_0000_map.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/q2DfcYPWBRRTX65/download","DCSM-FM_05nm_0001_map.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/TPeJjBPWDScKd8w/download","DCSM-FM_05nm_0002_map.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/pwKCWZ3jH2SRmgY/download","DCSM-FM_05nm_0003_map.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/p8H33Nzmpa3sREx/download","DCSM-FM_05nm_0004_map.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/cdSKmmgYgtwRQtA/download","DCSM-FM_05nm_0005_map.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/AmWenxpJSKJHNYQ/download","DCSM-FM_05nm_0006_map.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/xwgwo45pTmnD6cc/download","DCSM-FM_05nm_0007_map.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/cY47yHSot8rxmW5/download","DCSM-FM_05nm_0008_map.nc")
HTTP.download("https://nx7384.your-storageshare.de/s/y4se8798GsXfZEZ/download","DCSM-FM_05nm_0009_map.nc")

