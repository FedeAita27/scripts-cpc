!pip install cartopy
!pip install NetCDF4
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

windnetcdf4 = Dataset("C:/Users/Bolsista/Desktop/Federico/IC_CHICO/ERA5/Vetor vento (850 e 200 hPa)/era5_clim_february.nc")
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
clim_period_V = da_v.sel(valid_time=slice('1991-02-01', '2020-02-29')) # Sempre trocar de acordo com o mes
clim_period_U = da_u.sel(valid_time=slice('1991-02-01', '2020-02-29')) # Sempre trocar de acordo com o mes
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


anom_daily_U
anom_daily_V[1038]
anom_daily_U[1038]
anom_desired_V = anom_daily_V.sel(valid_time='2024-02-25')
anom_desired_U = anom_daily_U.sel(valid_time='2024-02-25')
anom_required_V = anom_desired_V.sel(pressure_level=850)
anom_required_U = anom_desired_U.sel(pressure_level=850)

wind_magnitude = np.sqrt(anom_required_U**2 + anom_required_V**2) # Cálculo de velocidade do vento



mask = wind_magnitude > 11
U_filtrado = anom_required_U.where(mask)
V_filtrado = anom_required_V.where(mask)

lon = wind_magnitude.longitude
lat = wind_magnitude.latitude

magnitude_levels = np.linspace(5,15,10)


fig, ax = plt.subplots(1, 1, figsize = (16, 8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
#ax.set_title('Anomalia de vetor vento em 850 hPa - 17/01/2024 (1991-2020)')

wind_contour = ax.contourf(lon ,lat, wind_magnitude, cmap='Blues',
                 levels=magnitude_levels,
                 extend='max', transform=ccrs.PlateCarree())

tiler = cimgt.GoogleTiles(style='satellite')
ax.add_image(tiler,6) # Aumentar o valor, aumenta a resolução

pular=12
ax.quiver(lon[::pular], lat[::pular], U_filtrado[::pular, ::pular].values,
          V_filtrado[::pular, ::pular].values,
           transform=ccrs.PlateCarree(), scale=400, color='black',
           headwidth=2, headlength=3, linewidth=8, alpha=1)
ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False),
              linewidth=1, color='gray', alpha=0.5, linestyle='--', transform=ccrs.PlateCarree())

ax.set_extent([-15, -105, 10, -90], crs=ccrs.PlateCarree())

# Barra de cores
colorbar_ticks = np.linspace(5,15,11)
colorbar = plt.colorbar(wind_contour, ticks=colorbar_ticks)
colorbar.set_label('m/s',size=20)
