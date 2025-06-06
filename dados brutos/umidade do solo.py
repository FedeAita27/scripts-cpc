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

nc = Dataset('C:/Users/feder/OneDrive/ERA5/netCDFs/gpcp_v02r03_monthly_d201402.nc')
ds = xr.open_dataset(xr.backends.NetCDF4DataStore(nc))
ds_180 = ds.assign_coords(longitude=(((ds.longitude + 180) % 360)-180)).sortby('longitude')
da = ds_180['precip']
da

plot_precip = da.sel(time='2014-02-01')

#fev = da_us.sel(valid_time=slice('2020-02-01', '2024-02-01'))
#mar = da_us.sel(valid_time=slice('2020-03-01','2024-03-01'))
#abr = da_us.sel(valid_time=slice('2020-04-01','2024-04-01'))
mai = da_us.sel(valid_time=slice('2020-05-01','2024-05-01'))
#media_dado = dado.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")

#media_fev = fev.groupby('valid_time.year').mean(keep_attrs=True).mean(dim='year')
#media_mar = mar.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")
#media_abr = abr.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")
media_mai = mai.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")
#media_dado = dado.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")

tittle_data = 'maio 2020-2024'

reg = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines',
                                   scale='50m', facecolor='none')

fig, ax = plt.subplots(1, 1, figsize = (16, 8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.add_feature(reg)
ax.set_extent([-60, -45, -35, -20])
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.set_title(f'Volume de água no solo - {tittle_data}', fontsize=16)
ax.gridlines(draw_labels=dict(left=False, bottom=False, top=False, right=False),
             linewidth=1, color='gray', alpha=0.5, linestyle='--')
magnitude_levels = np.linspace(0, 35, 13)

contour_temp = ax.contourf(ds_180.longitude, ds_180.latitude, plot_precip,
                           levels=magnitude_levels, extend='both', cmap='YlGnBu',
                           transform=ccrs.PlateCarree())

colorbar_ticks = np.linspace(0,35,7)
cbar = plt.colorbar(contour_temp, ax=ax, orientation='vertical', 
                    ticks=colorbar_ticks, pad=0.05, aspect=20)
cbar.set_label('mm', size=16)