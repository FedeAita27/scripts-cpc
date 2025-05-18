import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import scipy

windsnetcdf4 = ('C:/Users/feder/OneDrive/ERA5/TCC/nc/winds_may.nc')

ds_wind = xr.open_dataset(windsnetcdf4)
ds_wind
#ds_wind['date']
 
ds_180 = ds_wind.assign_coords(longitude=(((ds_wind.longitude + 180) % 360) - 180)).sortby('longitude')
da_v = ds_180['v']
da_u = ds_180['u']

da_v[0,:,:].plot()
def calculate_may_wind_anomalies(da_u, da_v, reference_period=('1991-05-01', '2020-05-01')):
    # 1. Filtrar apenas dados de setembro
    may_u = da_u.sel(valid_time=da_u['valid_time.month'] == 5)
    may_v = da_v.sel(valid_time=da_v['valid_time.month'] == 5)
    
    # 2. Calcular climatologia (média 1991-2020)
    clim_u = may_u.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    clim_v = may_v.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    
    # 3. Calcular anomalias (dados observados - climatologia)
    anom_u = may_u.groupby('valid_time.month') - clim_u
    anom_v = may_v.groupby('valid_time.month') - clim_v
 
    return {
        'anom_u': anom_u,
        'anom_v': anom_v,
        'clim_u': clim_u,
        'clim_v': clim_v
    }

# Calcular anomalias
results = calculate_may_wind_anomalies(da_u, da_v)

# Selecionar anos específicos em 200 e 850 hPa
target_date = '2024-05-01'
anom = results['anom_u'].sel(valid_time=target_date, pressure_level=200) 

fig, ax = plt.subplots(figsize=(16, 8), subplot_kw={'projection': ccrs.Orthographic(-50, -60)})

# Níveis e colormap
levels = np.linspace(-15, 15, 13)
cmap = plt.get_cmap('jet').copy()
cmap.set_over('darkred')
cmap.set_under('darkblue')

# Adicionar feições ao mapa
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)

# Contourf com foco na América do Sul
cf = ax.contourf(
    da_v.longitude, da_v.latitude, anom,
    levels=levels, cmap=cmap,
    extend='both', transform=ccrs.PlateCarree()
)

# Título e grade
year = pd.to_datetime(target_date).year
ax.set_title(f'Maio {year}', fontsize=16)
ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('m s⁻¹', fontsize=20)
cbar.ax.tick_params(labelsize=20)

plt.tight_layout()
plt.show()