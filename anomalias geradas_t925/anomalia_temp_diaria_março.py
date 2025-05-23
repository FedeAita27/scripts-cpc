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

tempnetcdf = Dataset('C:/Users/Bolsista/Desktop/Federico/IC_CHICO/ERA5/Temperatura (925 hPa)/era5_temp_clim_march.nc')
tempnetcdf

ds_temp = xr.open_dataset(xr.backends.NetCDF4DataStore(tempnetcdf))
ds_temp

da = ds_temp['t']
da

da[0,:,:].plot()

da_degc = da - 273.15
da_degc = da_degc.assign_attrs(da.attrs)
da_degc.attrs['units'] = '°C'

da_degc[0,:,:].plot()

yearly_mean= da.groupby('valid_time.year').mean(keep_attrs=True)
ref = yearly_mean.where((yearly_mean.year > 1990) & (yearly_mean.year < 2021), drop=True)
ref_mean = ref.mean(dim="year", keep_attrs=True)
clim_period = da.sel(valid_time=slice('1991-03-01', '2020-03-31')) # Sempre trocar de acordo com o mes
clim_month = clim_period.groupby('valid_time.month').mean()
weights= np.cos(np.deg2rad(clim_month.latitude))
weights.name = "weights"
clim_month_weighted = clim_month.weighted(weights)
mean = clim_month_weighted.mean(["longitude", "latitude"])
anom_daily = da.groupby('valid_time.month') - clim_month

anom_desired = anom_daily.sel(valid_time='2024-03-12', pressure_level=925)

lon = anom_desired.longitude
lat = anom_desired.latitude

magnitude_levels = np.linspace(-5, 5, 11)

#temp_magnitude = anom_desired > 

fig, ax = plt.subplots(1, 1, figsize = (10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.set_title('Anomalia de temperatura 925 hPa - 12/03/2024 (1991-2020)')

temp_contour = ax.contourf(lon ,lat, anom_desired, cmap='coolwarm',
                 levels=magnitude_levels,
                 extend='both', transform=ccrs.PlateCarree())

tiler = cimgt.GoogleTiles(style='satellite')
ax.add_image(tiler,6) # Aumentar o valor, aumenta a resolução



ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False),
              linewidth=1, color='gray', alpha=0.5, linestyle='--', transform=ccrs.PlateCarree())

#ax.set_extent([-20, -80, -50, 0], crs=ccrs.PlateCarree())   
ax.set_extent([-15, -105, 10, -90], crs=ccrs.PlateCarree())
ax.plot(-51.12, -30.02, marker='o', color='yellow', markersize=7, alpha=0.7)

# Barra de cores
colorbar_ticks = np.linspace(-5,5,11)
colorbar = plt.colorbar(temp_contour, ticks=colorbar_ticks)
colorbar.set_label('°C',size=20)
