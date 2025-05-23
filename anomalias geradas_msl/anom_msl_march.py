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
from scipy import ndimage
from netCDF4 import Dataset

mslnetcdf4 = Dataset('C:/Users/Bolsista/Desktop/Federico/IC_CHICO/ERA5/MSL/era5_msl_march.nc')

ds_msl = xr.open_dataset(xr.backends.NetCDF4DataStore(mslnetcdf4))
ds_msl

ds_msl_180 = ds_msl.assign_coords(longitude=(((ds_msl.longitude + 180) % 360) - 180)).sortby('longitude')

da_msl = ds_msl_180['msl']
da_msl
def Pa_to_mb(Pa):
    return Pa/100

msl_mb = Pa_to_mb(da_msl)

yearly_mean_msl = msl_mb.groupby('valid_time.year').mean(keep_attrs=True)
ref_winds_msl = yearly_mean_msl.where((yearly_mean_msl.year > 1990) & (yearly_mean_msl.year < 2021), drop=True)
ref_mean_winds_msl = ref_winds_msl.mean(dim="year", keep_attrs=True)
clim_period_msl = msl_mb.sel(valid_time=slice('1991-03-01', '2020-03-31')) # Sempre trocar de acordo com o mes
clim_month_msl = clim_period_msl.groupby('valid_time.month').mean()
weights_msl = np.cos(np.deg2rad(clim_month_msl.latitude))
weights_msl.name = "weights_msl"
clim_month_weighted_msl= clim_month_msl.weighted(weights_msl)
mean_msl = clim_month_weighted_msl.mean(["longitude", "latitude"])
anom_daily_msl = msl_mb.groupby('valid_time.month') - clim_month_msl

anom_desired = anom_daily_msl.sel(valid_time='2023-03-06')

lon = anom_desired.longitude
lat = anom_desired.latitude


min_msl = anom_desired.min().values
max_msl = anom_desired.max().values
magnitude_levels = np.linspace(-15, 15, 11)

reg = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines',
                                   scale='50m', facecolor='none')

fig, ax = plt.subplots(1, 1, figsize=(16,8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.add_feature(reg)
ax.set_title('Anomaly Sea level Pressure', size=12)
msl_contour = ax.contourf(lon,lat,anom_desired, cmap='Spectral_r',
            levels=magnitude_levels, linewidth=0.8, extend='both', transform=ccrs.PlateCarree())
tiler = cimgt.GoogleTiles(style='satellite')
ax.add_image(tiler,6) # Aumentar o valor, aumenta a resolução
ax.set_extent([-15, -105, 10, -80], crs=ccrs.PlateCarree())

#ax.clabel(msl_contour, inline=True, fontsize=8, fmt='%1.0f')

ax.plot(-50, -30, 'ro', markersize=8, transform=ccrs.PlateCarree())    
ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False),
    linewidth=1, color='gray', alpha=0.5, linestyle='--', transform = ccrs.PlateCarree())

colorbar_ticks = np.linspace(-15,15,11)
colorbar = plt.colorbar(msl_contour, ticks=colorbar_ticks)
colorbar.set_label('hPa',size=20)
