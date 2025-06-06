!pip install cartopy
!pip install netCDF4
!pip install scipy

import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import scipy

meannetcdf4 = Dataset('C:/Users/feder/Downloads/e6b2d6b1c37b0c63076d88e3ce4b786e.nc')
mean = xr.open_dataset(xr.backends.NetCDF4DataStore(meannetcdf4))

mean_180 = mean.assign_coords(longitude=(((mean.longitude + 180) % 360) - 180)).sortby('longitude')

msl = mean_180['msl']

def Pa_to_mb(Pa):
   return Pa / 100.0

msl_mb = Pa_to_mb(msl)

desired_date = '2024-05-05'
tittle_data = '05/05/2024'

plot_msl = msl_mb.sel(valid_time=desired_date)

fig, ax = plt.subplots(1, 1, figsize = (16, 8), subplot_kw={'projection': ccrs.Orthographic(-50, -30)})

#ax.stock_img()
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.set_title(f'Mean Sea Level Pressure - {tittle_data}', fontsize=16)
ax.gridlines(draw_labels=dict(left=False, bottom=False, top=False, right=False),
             linewidth=1, color='gray', alpha=0.5, linestyle='--')
magnitude_levels = np.linspace(970, 1030, 13)

# ISÓBARAS
#msl_contour = ax.contour(mean_180.longitude, mean_180.latitude, plot_msl,
 #   levels=magnitude_levels, extend='both', colors='k',
  #  transform=ccrs.PlateCarree())

# CONTORNO
msl_contour = ax.contourf(mean_180.longitude, mean_180.latitude, plot_msl,
    levels=magnitude_levels, extend='both', cmap='jet',
    transform=ccrs.PlateCarree())

#ax.clabel(msl_contour, inline=True, fontsize=8, fmt='%1.0f')

colorbar_ticks = np.linspace(970,1030,7)
cbar = plt.colorbar(msl_contour, ax=ax, orientation='vertical', 
                    ticks=colorbar_ticks, pad=0.05, aspect=20)
cbar.set_label('hPa', size=16)
