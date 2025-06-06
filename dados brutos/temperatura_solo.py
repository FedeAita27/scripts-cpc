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

t2mnc = Dataset('C:/Users/feder/OneDrive/ERA5/netCDFs/temp_solo_media_maio.nc')
t2m = xr.open_dataset(xr.backends.NetCDF4DataStore(t2mnc))
t2m_180 = t2m.assign_coords(longitude=(((t2m.longitude + 180) % 360)-180)).sortby('longitude')
da_t2m = t2m_180['stl1']

da_degc = da_t2m - 273.15
da_degc = da_degc.assign_attrs(da_t2m.attrs)
da_degc.attrs['units'] = '° C'
da_t2m = da_degc

#fev = da_t2m.sel(valid_time=slice('2014-02-01', '2024-02-29'))
#mar = da_t2m.sel(valid_time=slice('2020-03-01','2024-03-31'))
#abr = da_t2m.sel(valid_time=slice('2020-04-01','2024-04-30'))
mai = da_t2m.sel(valid_time=slice('2020-05-01','2024-05-31'))
#media_dado = dado.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")

#media_fev = fev.groupby('valid_time.year').mean(keep_attrs=True).mean(dim='year')
#media_mar = mar.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")
#media_abr = abr.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")
media_mai = mai.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")
#media_dado = dado.groupby('valid_time.year').mean(keep_attrs=True).mean(dim="year")

#desired_month = '2014-05-01'
tittle_data = 'maio 2020-2024'

#plot_t2m = da_t2m.sel(valid_time=desired_month)
reg = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines',
                                   scale='50m', facecolor='none')

fig, ax = plt.subplots(1, 1, figsize = (16, 8), subplot_kw={'projection': ccrs.PlateCarree()})


ax.add_feature(reg)
ax.set_extent([-60, -45, -35, -20])
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1)
ax.set_title(f'Temperatura média do solo {tittle_data}', fontsize=16)
ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False),
             linewidth=1, color='gray', alpha=0.5, linestyle='--')
magnitude_levels = np.linspace(0, 30, 13)

contour_temp = ax.contourf(t2m_180.longitude, t2m_180.latitude, media_mai,
                           levels=magnitude_levels, extend='both', cmap='coolwarm',
                           transform=ccrs.PlateCarree())

colorbar_ticks = np.linspace(0,30,7)
cbar = plt.colorbar(contour_temp, ax=ax, orientation='vertical', 
                    ticks=colorbar_ticks, pad=0.05, aspect=20)
cbar.set_label('°C', size=16)