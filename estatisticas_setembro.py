# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 16:36:13 2025

@author: feder
"""

import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import scipy
import cdsapi
import matplotlib.patches as mpatches
from matplotlib.path import Path



dataset = "reanalysis-era5-pressure-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "geopotential",
        "u_component_of_wind",
        "v_component_of_wind"
    ],
    "pressure_level": [
        "200", "250", "500",
        "850"
    ],
    "year": [
        "1982", "1983", "1984",
        "1985", "1986", "1987",
        "1988", "1989", "1990",       
        "1991", "1992", "1993",
        "1994", "1995", "1996",
        "1997", "1998", "1999",
        "2000", "2001", "2002",
        "2003", "2004", "2005",
        "2006", "2007", "2008",
        "2009", "2010", "2011",
        "2012", "2013", "2014",
        "2015", "2016", "2017",
        "2018", "2019", "2020",
        "2021", "2022", "2023"
    ],
    "month": ["09"],
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [-30, -180, -90, 180]
}

url='https://cds.climate.copernicus.eu/api' #TROCAR URL
key='a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a'  #TROCAR KEY

client = cdsapi.Client(url, key)
client.retrieve(dataset, request).download()

netcdf4 = Dataset('2d8a8c173de0968e1f44b3b856c54f31.nc')
ds = xr.open_dataset(xr.backends.NetCDF4DataStore(netcdf4))

ds_180 = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')

da_v = ds_180['v'].sel(pressure_level=850)

# climatologia de setembro de 1991 a 2020
climatologia_set = da_v.sel(valid_time=da_v['valid_time'].dt.year.isin(range(1991, 2021)))
climatologia= climatologia_set.mean(dim='valid_time')
climatologia_std = climatologia_set.std(dim='valid_time')

# Setembro de 2023
set_2023 = da_v.sel(valid_time=da_v['valid_time'].dt.year == 2023)

# Anomalia
anomalia_2023 = set_2023 - climatologia


# quão anômalo foi setembro de 2023 (z-score)
z_score = anomalia_2023 / climatologia_std
z_score_selected = z_score.sel(pressure_level=850)

# Desvio padrão
desv = anomalia_2023.std(dim='valid_time')

fig, ax = plt.subplots(figsize=(16, 8), subplot_kw={'projection': ccrs.Orthographic(-30, -90)})

ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)

cf = ax.contourf(da_v.longitude, da_v.latitude, anomalia_2023.squeeze(),
                 cmap="RdBu_r",
                 extend='both', transform=ccrs.PlateCarree())

contours = ax.contour(
    da_v.longitude,
    da_v.latitude,
    desv,
    levels=[-2, 2],  
    colors='black',
    linewidths=5,
    linestyles='--',
    transform=ccrs.PlateCarree()
)


ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='max')
cbar.set_label('Desvio padrão', fontsize=20)
cbar.ax.tick_params(labelsize=20)

plt.tight_layout()
plt.show()