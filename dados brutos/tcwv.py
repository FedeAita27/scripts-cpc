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


# TCWV 
tcwvnetcdf4 = Dataset('C:/Users/feder/OneDrive/ERA5/netCDFs/tcwv_dez_2024.nc')
ds_tcwv = xr.open_dataset(xr.backends.NetCDF4DataStore(tcwvnetcdf4))
ds_tcwv_180 = ds_tcwv.assign_coords(longitude=(((ds_tcwv.longitude + 180) % 360) - 180)).sortby('longitude')
da_tcwv = ds_tcwv_180['tcwv']

desired_time = '2024-12-26'
plot_tcwv = da_tcwv.sel(valid_time=desired_time)

# VECTOR WIND
uwindnc = Dataset('C:/Users/feder/OneDrive/ERA5/netCDFs/uwind_dez_2024.nc')
vwindnc = Dataset('C:/Users/feder/OneDrive/ERA5/netCDFs/vwind_dez_2024.nc')
uwind = xr.open_dataset(xr.backends.NetCDF4DataStore(uwindnc))
vwind = xr.open_dataset(xr.backends.NetCDF4DataStore(vwindnc))
uwind_180 = uwind.assign_coords(longitude=(((uwind.longitude + 180)% 360)-180)).sortby('longitude')
vwind_180 = vwind.assign_coords(longitude=(((vwind.longitude + 180)% 360)-180)).sortby('longitude')
da_uwind = uwind_180['u']
da_vwind = vwind_180['v']
'''
vectorwindnc = Dataset('C:/Users/feder/OneDrive/ERA5/netCDFs/vectorwinds_160124.nc')
vectorwind =  xr.open_dataset(xr.backends.NetCDF4DataStore(vectorwindnc))
vectorwind
vectorwind_180 = vectorwind.assign_coords(longitude=(((vectorwind.longitude + 180)% 360)-180)).sortby('longitude')
da_vwind = vectorwind_180['v']
da_uwind = vectorwind_180['u']
'''
uwind_850 = da_uwind.sel(valid_time=desired_time, pressure_level=850)
vwind_850 = da_vwind.sel(valid_time=desired_time, pressure_level=850)

wind_magnitude = np.sqrt(uwind_850**2 + vwind_850**2)

# MÁSCARA PARA PLOTAR OS VETORES APENAS ACIMA DO VALOR PEDIDO
mask = wind_magnitude >= 8
U_filtered = uwind_850.where(mask)
V_filtered = vwind_850.where(mask)


fig, ax = plt.subplots(1, 1, figsize=(20,10), 
          subplot_kw={'projection': ccrs.Orthographic(60, -70)})

ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.COASTLINE)
#ax.set_global([30, 110, -50, 20])  # Mantém o formato ortográfico sem recortar
ax.set_extent([20, 120, -70, 10], crs=ccrs.PlateCarree())
ax.gridlines(draw_labels=dict(
    left=False, bottom=False, right=False, top=False), 
    linewidth=1, color='gray', alpha=0.5, linestyle='--')
ax.set_title(f'Coluna total de vapor de água e vetor vento 850 hPa - {desired_time}', fontsize=16)
magnitude_levels = np.linspace(0, 50, 13)


contour_vapour = ax.contourf(ds_tcwv_180.longitude, ds_tcwv_180.latitude, plot_tcwv,
                             cmap='GnBu', extend='max', levels=magnitude_levels, transform=ccrs.PlateCarree())

pular = 10
ax.quiver(uwind_850.longitude.values[::pular], uwind_850.latitude.values[::pular],
          U_filtered.values[::pular, ::pular], V_filtered.values[::pular, ::pular],
          color='black', scale=500, headwidth=2, headlength=3, linewidth=1, alpha=1,
          transform=ccrs.PlateCarree())

colorbar_ticks = np.linspace(0, 50, 7) 
colorbar = plt.colorbar(contour_vapour, orientation='vertical', ticks=colorbar_ticks) 
colorbar.set_label('mm',size=20)