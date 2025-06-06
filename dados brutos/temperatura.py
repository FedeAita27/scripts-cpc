!pip install cartopy
!pip install netCDF4

import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import pandas as pd
from netCDF4 import Dataset

t925nc = Dataset('C:/Users/feder/OneDrive/ERA5/netCDFs/3379107b325a645f91f1a06d612f0f24.nc')
t925 = xr.open_dataset(xr.backends.NetCDF4DataStore(t925nc))
t925_180 = t925.assign_coords(longitude=(((t925.longitude + 180) % 360)-180)).sortby('longitude')

da_t925 = t925_180['t']

da_degc = da_t925 - 273.15
da_degc = da_degc.assign_attrs(da_t925.attrs)
da_degc.attrs['units'] = '° C'
da_t925 = da_degc

desired_data = '2024-05-05'
tittle_data = '05/05/2024'

plot_t925 = da_t925.sel(valid_time=desired_data, pressure_level=925)

fig, ax = plt.subplots(1, 1, figsize = (16, 8), subplot_kw={'projection': ccrs.Orthographic(-50, -30)})

ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.set_title(f'Temperatura a 925 hPa - {tittle_data}', fontsize=16)
ax.gridlines(draw_labels=dict(left=False, bottom=False, top=False, right=False),
             linewidth=1, color='gray', alpha=0.5, linestyle='--')
magnitude_levels = np.linspace(-30, 30, 13)

contour_temp = ax.contourf(t925_180.longitude, t925.latitude, plot_t925,
                           levels=magnitude_levels, extend='both', cmap='coolwarm',
                           transform=ccrs.PlateCarree())

colorbar_ticks = np.linspace(-30,30,7)
cbar = plt.colorbar(contour_temp, ax=ax, orientation='vertical', 
                    ticks=colorbar_ticks, pad=0.05, aspect=20)
cbar.set_label('°C', size=16)