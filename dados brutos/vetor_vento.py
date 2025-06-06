!pip install cartopy
!pip install netCDF4

import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# DAILY MEAN
uwindnc = Dataset('C:/Users/feder/OneDrive/ERA5/u_wind_3103.nc')
vwindnc = Dataset('C:/Users/feder/OneDrive/ERA5/v_wind_3103.nc')
uwind = xr.open_dataset(xr.backends.NetCDF4DataStore(uwindnc))
vwind = xr.open_dataset(xr.backends.NetCDF4DataStore(vwindnc))
uwind_180 = uwind.assign_coords(longitude=(((uwind.longitude + 180)% 360)-180)).sortby('longitude')
vwind_180 = vwind.assign_coords(longitude=(((vwind.longitude + 180)% 360)-180)).sortby('longitude')
da_uwind = uwind_180['u']
da_vwind = vwind_180['v']



# HOURLY DATA
vectorwindnc = Dataset('C:/Users/feder/OneDrive/ERA5/netCDFs/vw_3103.nc')
vectorwind =  xr.open_dataset(xr.backends.NetCDF4DataStore(vectorwindnc))
vectorwind
vectorwind_180 = vectorwind.assign_coords(longitude=(((vectorwind.longitude + 180)% 360)-180)).sortby('longitude')
da_vwind = vectorwind_180['v']
da_uwind = vectorwind_180['u']

desired_data = '2025-03-31'

# VETOR VENTO EM 850 hPa
uwind_850 = da_uwind.sel(valid_time=desired_data, pressure_level=850)
vwind_850 = da_vwind.sel(valid_time=desired_data, pressure_level=850)

# VETOR VENTO EM 200 hPa
uwind_200 = da_uwind.sel(valid_time=desired_data, pressure_level=200)
vwind_200 = da_vwind.sel(valid_time=desired_data, pressure_level=200)

wind_magnitude = np.sqrt(uwind_850**2 + vwind_850**2)

mask = wind_magnitude > 5
U_filtered = uwind_850.where(mask)
V_filtered = vwind_850.where(mask)

fig, ax = plt.subplots(1, 1, figsize = (16, 8), subplot_kw={'projection': ccrs.Orthographic(-50, -50)})

#ax.stock_img()
#ax.set_extent([-30, -80, -5, -45])
#ax.plot(-51.23, -30.03, '^', color='black' ,markersize=10)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.set_title(f'Vetor vento em 850 hPa - {desired_data}', fontsize=16)
ax.gridlines(draw_labels=dict(left=False, bottom=False, top=False, right=False),
             linewidth=1, color='gray', alpha=0.5, linestyle='--')
magnitude_levels = np.linspace(4, 24, 13)

contour_wind = ax.contourf(uwind_180.longitude, uwind_180.latitude, wind_magnitude,
                           levels=magnitude_levels, extend='max', cmap='GnBu',
                           transform=ccrs.PlateCarree())
pular = 8
ax.quiver(uwind_850.longitude.values[::pular], uwind_850.latitude.values[::pular],
          U_filtered.values[::pular, ::pular], V_filtered.values[::pular, ::pular],
          color='black', scale=600, headwidth=2, headlength=3, linewidth=1, alpha=1,
          transform=ccrs.PlateCarree())

colorbar_ticks = np.linspace(4,24,7)
cbar = plt.colorbar(contour_wind, orientation='vertical', ticks=colorbar_ticks)
cbar.set_label('m/s', size=16)
