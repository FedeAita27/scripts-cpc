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
import cdsapi
import scipy
from scipy.ndimage import minimum_filter, maximum_filter, label, find_objects


dataset = "derived-era5-single-levels-daily-statistics"
request = {
    "product_type": "reanalysis",
    "variable": ["mean_sea_level_pressure"],
    "year": "2023",
    "month": ["09"],
    "day": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12",
        "13", "14", "15"
    ],
    "daily_statistic": "daily_mean",
    "time_zone": "utc+03:00",
    "frequency": "1_hourly"
}

url='https://cds.climate.copernicus.eu/api'
key='a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a'

client = cdsapi.Client(url,key)
client.retrieve(dataset, request).download()

meannetcdf4 = Dataset('e414c2d970f74ac71e4eff188501f509.nc')
mean = xr.open_dataset(xr.backends.NetCDF4DataStore(meannetcdf4))

mean_180 = mean.assign_coords(longitude=(((mean.longitude + 180) % 360) - 180)).sortby('longitude')

msl = mean_180['msl']

def Pa_to_mb(Pa):
   return Pa / 100.0

def plot_pressure_centers(pressure, lons, lats, ax, extent ,size=20):
    min_filt = minimum_filter(pressure, size=size, mode='reflect')
    minima = (pressure == min_filt)
    max_filt = maximum_filter(pressure, size=size, mode='reflect')
    maxima = (pressure == max_filt)

    # Baixas pressões (L)
    labeled, _ = label(minima)
    slices = find_objects(labeled)
    for dyx in slices:
        y, x = dyx[0].start, dyx[1].start
        ax.text(lons[y, x], lats[y, x], 'L', color='red', fontsize=14,
                fontweight='bold', ha='center', va='center', transform=ccrs.PlateCarree())

    # Altas pressões (H)
    labeled, _ = label(maxima)
    slices = find_objects(labeled)
    for dyx in slices:
        y, x = dyx[0].start, dyx[1].start
        ax.text(lons[y, x], lats[y, x], 'H', color='blue', fontsize=14,
                fontweight='bold', ha='center', va='center', transform=ccrs.PlateCarree())

msl_mb = Pa_to_mb(msl)


desired_date = '2023-09-14'

plot_msl = msl_mb.sel(valid_time=desired_date)

fig, ax = plt.subplots(1, 1, figsize = (10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.stock_img()
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.plot(51.3, 30.2, 'o', color='red')
#ax.set_title(f'Mean Sea Level Pressure - {tittle_data}', fontsize=16)
ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False),
             linewidth=1, color='gray', alpha=0.5, linestyle='--')
magnitude_levels = np.linspace(970, 1030, 20)
map_extent = [-80, -15, -80, 0]
ax.set_extent(map_extent, crs=ccrs.PlateCarree())


# ISÓBARAS
msl_contour = ax.contour(mean_180.longitude, mean_180.latitude, plot_msl,
    levels=magnitude_levels, extend='both', colors='k',
    transform=ccrs.PlateCarree())

# CONTORNO
#msl_contour = ax.contourf(mean_180.longitude, mean_180.latitude, plot_msl,
 #   levels=magnitude_levels, extend='both', cmap='jet',
  #  transform=ccrs.PlateCarree())

ax.clabel(msl_contour, inline=True, fontsize=8, fmt='%1.0f')

'''
colorbar_ticks = np.linspace(970,1030,7)
cbar = plt.colorbar(msl_contour, ax=ax, orientation='vertical', 
                    ticks=colorbar_ticks, pad=0.05, aspect=20)
cbar.set_label('hPa', size=16)
'''