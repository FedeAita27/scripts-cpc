from netCDF4 import Dataset 
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

gptnetcdf4 = Dataset('C:/Users/Bolsista/Desktop/Federico/IC_CHICO/ERA5/Anomalias_diarias/Geopotencial/era_gpt_march.nc')

ds_gpt = xr.open_dataset(xr.backends.NetCDF4DataStore(gptnetcdf4))
ds_gpt_180 = ds_gpt.assign_coords(longitude=(((ds_gpt.longitude + 180) % 360) - 180)).sortby('longitude')

da_gpt = ds_gpt_180['z']

def gpt_conv(z):
    return z/9.80665

gpt = gpt_conv(da_gpt)

yearly_mean_gpt = gpt.groupby('valid_time.year').mean(keep_attrs=True)
ref_winds_gpt = yearly_mean_gpt.where((yearly_mean_gpt.year > 1990) & (yearly_mean_gpt.year < 2021), drop=True)
ref_mean_winds_gpt = ref_winds_gpt.mean(dim="year", keep_attrs=True)
clim_period_gpt = gpt.sel(valid_time=slice('1991-03-01', '2020-03-31')) # Sempre trocar de acordo com o mes
clim_month_gpt = clim_period_gpt.groupby('valid_time.month').mean()
weights_gpt = np.cos(np.deg2rad(clim_month_gpt.latitude))
weights_gpt.name = "weights_gpt"
clim_month_weighted_gpt= clim_month_gpt.weighted(weights_gpt)
mean_gpt = clim_month_weighted_gpt.mean(["longitude", "latitude"])
anom_daily_gpt = gpt.groupby('valid_time.month') - clim_month_gpt

anom_desired = anom_daily_gpt.sel(valid_time='2024-03-16', pressure_level=925)

lon = anom_desired.longitude
lat = anom_desired.latitude

min_gpt = anom_desired.min().values
max_gpt = anom_desired.max().values

plot_level = 925
data = '16/03/2024'
fig, ax = plt.subplots(1, 1, figsize=(16,8), subplot_kw={'projection': ccrs.Orthographic(-50, -80)})

ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
#ax.add_feature(reg)
ax.set_title(f'Anomaly geopotential {plot_level} hPa - {data}', size=12)
gpt_contour = ax.contourf(lon,lat,anom_desired, cmap='jet',
            levels=np.linspace(-200, 200, 17), linewidth=0.8, extend='both', transform=ccrs.PlateCarree())
tiler = cimgt.GoogleTiles(style='satellite')
ax.add_image(tiler,6) # Aumentar o valor, aumenta a resolução

#ax.clabel(msl_contour, inline=True, fontsize=8, fmt='%1.0f')

#ax.plot(-50, -30, 'ro', markersize=8, transform=ccrs.PlateCarree())    
ax.gridlines(draw_labels=False,linewidth=1, color='gray', 
             alpha=0.5, linestyle='--', transform = ccrs.Orthographic())

colorbar_ticks = np.linspace(-200,200,9)
colorbar = plt.colorbar(gpt_contour, ticks=colorbar_ticks)
colorbar.set_label('m',size=20)
