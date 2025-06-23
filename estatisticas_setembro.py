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

url='https://cds.climate.copernicus.eu/api'
key='a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a'

client = cdsapi.Client(url, key)
client.retrieve(dataset, request).download()

netcdf4 = Dataset('e372092be0038c9869bbd90fa30ba79c.nc')
ds = xr.open_dataset(xr.backends.NetCDF4DataStore(netcdf4))

ds_180 = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby('longitude')

da_v = ds_180['v']

# Climatologia 
climatologia = da_v.sel(valid_time=slice('1991-09-01', '2020-09-01')).mean('valid_time')

# Anos específicos
anos = da_v.sel(valid_time=da_v['valid_time'].dt.year.isin([1982, 1997, 2015, 2023]))

# Anomalias
anomalias = anos - climatologia

# Desvio padrão
desv = anomalias.std(dim='valid_time')
desv_850 = desv.sel(pressure_level=850)


desv_selected = desv.sel(pressure_level=850)

vmin = -3
vmax = 3

levels = np.linspace(vmin, vmax, 13)

fig, ax = plt.subplots(figsize=(16, 8), subplot_kw={'projection': ccrs.Orthographic(-30, -90)})

cmap = plt.get_cmap('jet').copy()

ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)

cf = ax.contourf(da_v.longitude, da_v.latitude, desv_selected,
                 levels=levels, cmap=cmap,
                 extend='both', transform=ccrs.PlateCarree())

ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='max')
cbar.set_label('Desvio padrão', fontsize=20)
cbar.ax.tick_params(labelsize=20)

plt.tight_layout()
plt.show()