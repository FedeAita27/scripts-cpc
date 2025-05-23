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

windnetcdf4 = Dataset("C:/Users/Aluno/Desktop/Federico/Graduacao_Federico/IC_CHICO/ERA5/Anomalias_diarias/Vetor vento (850 e 200 hPa)/era5_clim_dezembro.nc")
windnetcdf4

ds_wind=xr.open_dataset(xr.backends.NetCDF4DataStore(windnetcdf4))
ds_wind

#ds_wind_180 = ds_wind.assign_coords(longitude=(((ds_wind.longitude + 180) % 360) - 180)).sortby('longitude')


da_v = ds_wind['v']
da_u = ds_wind['u']


da_v[0,:,:].plot()
da_u[0,:,:].plot()


yearly_mean_V = da_v.groupby('valid_time.year').mean(keep_attrs=True)
yearly_mean_U = da_u.groupby('valid_time.year').mean(keep_attrs=True)
ref_winds_V = yearly_mean_V.where((yearly_mean_V.year > 1990) & (yearly_mean_V.year < 2021), drop=True)
ref_winds_U = yearly_mean_U.where((yearly_mean_U.year > 1990) & (yearly_mean_U.year < 2021), drop=True)
ref_mean_winds_V = ref_winds_V.mean(dim="year", keep_attrs=True)
ref_mean_winds_U = ref_winds_U.mean(dim="year", keep_attrs=True)
clim_period_V = da_v.sel(valid_time=slice('1991-12-01', '2020-12-31')) # Sempre trocar de acordo com o mes
clim_period_U = da_u.sel(valid_time=slice('1991-12-01', '2020-12-31')) # Sempre trocar de acordo com o mes
clim_month_V = clim_period_V.groupby('valid_time.month').mean()
clim_month_U = clim_period_U.groupby('valid_time.month').mean()
weights_V = np.cos(np.deg2rad(clim_month_V.latitude))
weights_U = np.cos(np.deg2rad(clim_month_U.latitude))
weights_V.name = "weights_V"
weights_U.name = "weights_U"
clim_month_weighted_V = clim_month_V.weighted(weights_V)
clim_month_weighted_U = clim_month_U.weighted(weights_U)
mean_V = clim_month_weighted_V.mean(["longitude", "latitude"])
mean_U = clim_month_weighted_U.mean(["longitude", "latitude"])
anom_daily_V = da_v.groupby('valid_time.month') - clim_month_V
anom_daily_U = da_u.groupby('valid_time.month') - clim_month_U


desired_data = '2024-12-24'

anom_desired_V = anom_daily_V.sel(valid_time=desired_data)
anom_desired_U = anom_daily_U.sel(valid_time=desired_data)
anom_required_V = anom_desired_V.sel(pressure_level=850)
anom_required_U = anom_desired_U.sel(pressure_level=850)


wind_magnitude = np.sqrt(anom_required_U**2 + anom_required_V**2) # Cálculo de velocidade do vento


mask = wind_magnitude > 5
U_filtrado = anom_required_U.where(mask)
V_filtrado = anom_required_V.where(mask)

lon = wind_magnitude.longitude
lat = wind_magnitude.latitude

magnitude_levels = np.linspace(5,35,13)

fig, ax = plt.subplots(1, 1, figsize = (16, 8), subplot_kw={'projection': ccrs.Orthographic(60, -70)})

ax.stock_img()
ax.set_extent([20, 120, -70, 10], crs=ccrs.PlateCarree())
#ax.plot(100.42, -65.31, "^", color='red', markersize=10, transform=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.set_title(f'Anomalia de vetor vento em 850 hPa - {desired_data} ')

wind_contour = ax.contourf(lon ,lat, wind_magnitude, cmap='jet',
                 levels=magnitude_levels,
                 extend='max', transform=ccrs.PlateCarree())


pular=10
ax.quiver(lon[::pular], lat[::pular], U_filtrado[::pular, ::pular].values,
          V_filtrado[::pular, ::pular].values,
           transform=ccrs.PlateCarree(), scale=500, color='black',
           headwidth=2, headlength=3, linewidth=8, alpha=1)
ax.gridlines(draw_labels=dict(left=False, bottom=False, top=False, right=False),
              linewidth=1, color='gray', alpha=0.5, linestyle='--', transform=ccrs.PlateCarree())

#ax.set_extent([-20, -80, -50, 0], crs=ccrs.PlateCarree())   
#ax.set_extent([-15, -105, 10, -90], crs=ccrs.PlateCarree())
#ax.plot(-51.12, -30.02, marker='o', color='black', markersize=10, alpha=0.7)

# Barra de cores
colorbar_ticks = np.linspace(5,35,7)
colorbar = plt.colorbar(wind_contour, ticks=colorbar_ticks)
colorbar.set_label('m/s',size=20)
colorbar.ax.set_yticklabels([f'{x:.0f}' for x in colorbar_ticks])

