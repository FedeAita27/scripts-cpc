!pip install cartopy
!pip install matplotlib
!pip install netCDF4
!pip install xarray

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

windnetcdf4 = Dataset("C:/Users/Federico/Downloads/60449d9760b70d984871fa2f7e6f6a1b.nc")
windnetcdf4


ds_wind=xr.open_dataset(xr.backends.NetCDF4DataStore(windnetcdf4))
ds_wind

ds_wind_180 = ds_wind.assign_coords(longitude=(((ds_wind.longitude + 180) % 360) - 180)).sortby('longitude')


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
clim_period_V = da_v.sel(valid_time=slice('1991-09-01', '2020-09-30')) # Sempre trocar de acordo com o mes
clim_period_U = da_u.sel(valid_time=slice('1991-09-01', '2020-09-30')) # Sempre trocar de acordo com o mes
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
anom_daily_V = da_v.groupby('valid_time.month') - mean_V
anom_daily_U = da_u.groupby('valid_time.month') - mean_U


anom_daily_U
anom_daily_V[1038]
anom_daily_U[1038]

anom_desired_V = anom_daily_V.sel(valid_time='2023-09-14', pressure_level=200)
anom_desired_U = anom_daily_U.sel(valid_time='2023-09-14', pressure_level=200)


#title_data = '2023-09-03'

wind_magnitude = np.sqrt(anom_desired_U**2 + anom_desired_V**2) # Cálculo de velocidade do vento


mask = wind_magnitude > 10
U_filtrado = anom_desired_U.where(mask)
V_filtrado = anom_desired_V.where(mask)

lon = wind_magnitude_desired.longitude
lat = wind_magnitude_desired.latitude


magnitude_levels = np.linspace(10,70,13)


fig, ax = plt.subplots(1,1,figsize = (10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
#ax.set_title(f"Anomalia de vetor vento - {title_data}")

wind_contour = ax.contourf(lon ,lat, wind_magnitude_desired, cmap='jet',
                 levels=magnitude_levels,
                 extend='both', transform=ccrs.PlateCarree())
#tiler = cimgt.GoogleTiles(style='satellite')
#ax.add_image(tiler,6) # Aumentar o valor, aumenta a resolução

pular=8
ax.quiver(lon[::pular], lat[::pular], U_filtrado[::pular, ::pular].values,
          V_filtrado[::pular, ::pular].values,
           transform=ccrs.PlateCarree(), scale=1000, color='black',
           headwidth=2, headlength=3, linewidth=8, alpha=1)
ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False),
              linewidth=1, color='gray', alpha=0.5, linestyle='--', transform=ccrs.PlateCarree())

ax.set_extent([-20, -80, -50, 0], crs=ccrs.PlateCarree())   
#ax.set_extent([-15, -105, 10, -90], crs=ccrs.PlateCarree())
ax.plot(-51.12, -30.02, marker='o', color='red', markersize=7, alpha=0.7)

# Barra de cores
colorbar_ticks = np.linspace(10,70,7)
colorbar = plt.colorbar(wind_contour, ticks=colorbar_ticks, orientation='horizontal', pad=0.1, aspect=40, shrink=0.7)
colorbar.set_label('m s-¹',size=18)

plt.show()