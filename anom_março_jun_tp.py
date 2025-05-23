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

tpnetcdf4= Dataset('C:/Users/Bolsista/Desktop/Federico/IC_CHICO/ERA5/Anomalias_diarias/Total precipitation/era5_tp_march.nc')
tpnetcdf4

ds_tp = xr.open_dataset(xr.backends.NetCDF4DataStore(tpnetcdf4))
ds_tp

ds_tp_180 = ds_tp.assign_coords(longitude=(((ds_tp.longitude +180) % 360) - 180)).sortby('longitude')
ds_tp_180['valid_time']

tp_mm = ds_tp_180['tp'] / 1000

tp_mm[0,:,:].plot()


yearly_mean_tp = tp_mm.groupby('valid_time.year').mean(keep_attrs=True)
ref_winds_tp = yearly_mean_tp.where((yearly_mean_tp.year > 1990) & (yearly_mean_tp.year < 2021), drop=True)
ref_mean_winds_tp = ref_winds_tp.mean(dim="year", keep_attrs=True)
clim_period_tp = tp_mm.sel(valid_time=slice('1991-03-01', '2020-03-31')) # Sempre trocar de acordo com o mes
clim_month_tp = clim_period_tp.groupby('valid_time.month').mean()
weights_tp = np.cos(np.deg2rad(clim_month_tp.latitude))
weights_tp.name = "weights_tp"
clim_month_weighted_tp= clim_month_tp.weighted(weights_tp)
mean_tp = clim_month_weighted_tp.mean(["longitude", "latitude"])
anom_daily_tp = tp_mm.groupby('valid_time.month') - clim_month_tp

anom_desired = anom_daily_tp.sel(valid_time='2023-03-06')

lon = anom_desired.longitude
lat = anom_desired.latitude

tp_mm = tp_mm.assign_coords(tp_mm.attrs)
tp_mm.attrs['units'] = 'mm'

magnitude_levels = np.linspace(0, 1, 11)

reg = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines',
                                   scale='50m', facecolor='none')

fig, ax = plt.subplots(1, 1, figsize=(16,8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.add_feature(reg)
ax.set_title('Anomaly total precipitation', size=12)
msl_contour = ax.contourf(lon,lat,anom_desired, cmap='jet',
            levels=magnitude_levels, linewidth=0.8, extend='max', transform=ccrs.PlateCarree())
tiler = cimgt.GoogleTiles(style='satellite')
ax.add_image(tiler,6) # Aumentar o valor, aumenta a resolução
ax.set_extent([-15, -105, 10, -80], crs=ccrs.PlateCarree())

#ax.clabel(msl_contour, inline=True, fontsize=8, fmt='%1.0f')

ax.plot(-50, -30, 'ro', markersize=8, transform=ccrs.PlateCarree())    
ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False),
    linewidth=1, color='gray', alpha=0.5, linestyle='--', transform = ccrs.PlateCarree())

colorbar_ticks = np.linspace(0,1,11)
colorbar = plt.colorbar(msl_contour, ticks=colorbar_ticks)
colorbar.set_label('mm',size=20)


