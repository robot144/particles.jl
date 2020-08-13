#! /bin/bash
#
# The swan data is not kept in the repository because of its size. The original data was downloaded with:
# wget -O swan_20171007_28.nc 'http://matroos.deltares.nl:80//matroos/scripts/matroos.pl?source=swan_dcsm&anal=000000000000&z=0&xmin=3.000000&xmax=10.000000&ymin=50.000000&ymax=58.000000&coords=WGS84&xmin_abs=3&xmax_abs=10&ymin_abs=50&ymax_abs=58&color=wave_dir_th0,wave_height_hm0,wave_period_tm10&interpolate=&now=202001090000&to=201710280000&from=201710070000&outputformat=nc&stridex=&stridey=&stridetime=1&xn=&yn=&celly=&cellx=&fieldoutput=wave_dir_th0,wave_height_hm0,wave_period_tm10&format=nc'

# The data is now cached in a different location to make sure it stays available in identical form.
# 

# 170Mb wave data from SWAN model
wget -O swan_20171007_28.nc  https://nx7384.your-storageshare.de/s/WPcFEfgtqA5bnsS/download

