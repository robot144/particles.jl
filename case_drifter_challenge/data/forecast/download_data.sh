#! /bin/bash
# This script downloads the relevant data for this case. This data is not inside the repository 
# because of its size.

# CMEMS flow field
# 84 Mb per file 
wget -O u0_2021-10-29_00-00-00_2021-11-07_00-00-00.nc https://nx7384.your-storageshare.de/s/kiq3tqKAWsRp8zb/download & 
wget -O v0_2021-10-29_00-00-00_2021-11-07_00-00-00.nc https://nx7384.your-storageshare.de/s/QpyqXfrj3P5rMSG/download

