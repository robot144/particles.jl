#!/usr/bin/env python
# download ERA5 surface data using the CDS API.

import cdsapi

c = cdsapi.Client()

# area:  starting latitude and longitude followed by the ending latitude and longitude
c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type':'reanalysis',
        'format':'netcdf',
        'variable':[
            '10m_u_component_of_wind','10m_v_component_of_wind','mean_sea_level_pressure',
            'mean_direction_of_wind_waves','mean_wave_direction', 'mean_wave_period',
            'peak_wave_period', 'significant_height_of_combined_wind_waves_and_swell', 
            'significant_height_of_wind_waves', 'u_component_stokes_drift', 'v_component_stokes_drift', 
        ],
        'area':'57/5/50/10',
        'year':'2017',
        'month':['03','04'],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09','10',
            '11','12','13',
            '14','15','16',
            '17','18','19','20'
            '21','22','23',
            '24','25','26',
            '27','28','29','20','31'
        ],
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ]
    },
    'era5_wind_201703_04.nc')
