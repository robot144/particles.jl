#! /bin/bash
# This script downloads the relevant data for this case. This data is not inside the repository 
# because of its size.

# CMEMS flow field
# 84 Mb per file 
wget -O u0_2021-10-29_00-00-00_2021-11-07_00-00-00.nc https://nx7384.your-storageshare.de/s/kiq3tqKAWsRp8zb/download & 
wget -O v0_2021-10-29_00-00-00_2021-11-07_00-00-00.nc https://nx7384.your-storageshare.de/s/QpyqXfrj3P5rMSG/download

# NOAA GFS wind fields - forecast mode
# 26.8 Mb per file
wget -O GFS_FFT_0.25_fc_202111060000.nc https://nx7384.your-storageshare.de/s/DJZd6q5cEqcPMf7/download
wget -O GFS_FFT_0.25_fc_202111070000.nc https://nx7384.your-storageshare.de/s/8ibMiAD24aQozY2/download
wget -O GFS_FFT_0.25_fc_202111080000.nc https://nx7384.your-storageshare.de/s/7cgo4LFCfeztRjT/download
wget -O GFS_FFT_0.25_fc_202111090000.nc https://nx7384.your-storageshare.de/s/BKDparkFkaLe9cb/download
wget -O GFS_FFT_0.25_fc_202111100000.nc https://nx7384.your-storageshare.de/s/BkojicYLW9QGALA/download
wget -O GFS_FFT_0.25_fc_202111110000.nc https://nx7384.your-storageshare.de/s/XoRRjdxoP6TBbY6/download
wget -O GFS_FFT_0.25_fc_202111120000.nc https://nx7384.your-storageshare.de/s/j2xsMX7L7KcP5Sz/download
wget -O GFS_FFT_0.25_fc_202111130000.nc https://nx7384.your-storageshare.de/s/WcdczKoeg7Z39gW/download
wget -O GFS_FFT_0.25_fc_202111140000.nc https://nx7384.your-storageshare.de/s/6yAyxkr25RCeR8i/download
wget -O GFS_FFT_0.25_fc_202111150000.nc https://nx7384.your-storageshare.de/s/NfWCiE8iWitxGKH/download
wget -O GFS_FFT_0.25_fc_202111160000.nc https://nx7384.your-storageshare.de/s/iY4gLpBt2pikFf7/download
wget -O GFS_FFT_0.25_fc_202111170000.nc https://nx7384.your-storageshare.de/s/YNpXPL9xp66e8HH/download

