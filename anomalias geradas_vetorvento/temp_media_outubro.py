!pip install cartopy
!pip install NetCDF4
!pip install scipy

import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import xarray as xr
import numpy as np
import pandas as pd
import scipy
from netCDF4 import Dataset

tempnetcdf4 = ('C:/Users/Bolsista/Desktop/Federico/IC_CHICO/nc/media_mensal_outubro.nc')

ds_temp = xr.open_dataset(tempnetcdf4)

ds_180 = ds_temp.assign_coords(longitude=(((ds_temp.longitude + 180) % 360) - 180)).sortby('longitude')
da_temp = ds_180['t2m']

da_degc = da_temp - 273.15
da_degc = da_degc.assign_attrs(da_temp.attrs)
da_degc.attrs['units'] = '°C'

yearly_mean= da_degc.groupby('valid_time.year').mean(keep_attrs=True)
ref = yearly_mean.where((yearly_mean.year > 1990) & (yearly_mean.year < 2021), drop=True)
ref_mean = ref.mean(dim="year", keep_attrs=True)
clim_period = da_degc.sel(valid_time=slice('1991-10-01', '2020-10-31')) # Sempre trocar de acordo com o mes
clim_month = clim_period.groupby('valid_time.month').mean()

clim_month_out = clim_month.sel(month=10)
clim_month

fig, ax = plt.subplots(1, 1, figsize = (10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

lon = clim_month.longitude
lat = clim_month.latitude

magnitude_levels = np.linspace(8, 30, 13)

ax.set_extent([-74, -34, -35, 5], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.add_feature(cfeature.STATES, linewidth=1)
#ax.set_title('temperatura méda outubro - (1991-2020)')

temp_contour = ax.contourf(lon, lat, clim_month_out, cmap='coolwarm',
                 levels=magnitude_levels,
                 extend='both', transform=ccrs.PlateCarree())


ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False),
              linewidth=1, color='gray', alpha=0.5, linestyle='--', transform = ccrs.PlateCarree())


# Barra de cores
colorbar_ticks = np.linspace(8,30,7)
colorbar = plt.colorbar(temp_contour, ticks=colorbar_ticks)
colorbar.set_label('°C',size=20)
