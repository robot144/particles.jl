#! /bin/bash

wget -O swan_20200110_11.nc 'https://matroos.rws.nl:80//direct/get_matroos.php?source=swan_dcsm&anal=000000000000&z=0&xmin=-12.000000&xmax=9.000000&ymin=48.000000&ymax=64.099998&coords=WGS84&xmin_abs=-12&xmax_abs=9&ymin_abs=48&ymax_abs=64.099998474121&color=wave_dir_th0,wave_height_hm0,wave_period_tm10&interpolate=&now=202001090000&to=202001110000&from=202001100000&outputformat=nc&stridex=&stridey=&stridetime=1&xn=&yn=&celly=&cellx=&fieldoutput=wave_dir_th0,wave_height_hm0,wave_period_tm10&format=nc'
