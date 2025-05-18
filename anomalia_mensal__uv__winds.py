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

windsnetcdf4 = ('C:/Users/feder/OneDrive/ERA5/TCC/nc/winds_september.nc')

ds_wind = xr.open_dataset(windsnetcdf4)
ds_wind
#ds_wind['date']
 
ds_180 = ds_wind.assign_coords(longitude=(((ds_wind.longitude + 180) % 360) - 180)).sortby('longitude')
da_v = ds_180['v']
da_u = ds_180['u']

da_v[0,:,:].plot()
def calculate_september_wind_anomalies(da_u, da_v, reference_period=('1991-09-01', '2020-09-01')):
    # 1. Filtrar apenas dados de setembro
    sept_u = da_u.sel(valid_time=da_u['valid_time.month'] == 9)
    sept_v = da_v.sel(valid_time=da_v['valid_time.month'] == 9)
    
    # 2. Calcular climatologia (média 1991-2020)
    clim_u = sept_u.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    clim_v = sept_v.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    
    # 3. Calcular anomalias (dados observados - climatologia)
    anom_u = sept_u.groupby('valid_time.month') - clim_u
    anom_v = sept_v.groupby('valid_time.month') - clim_v
 
    return {
        'anom_u': anom_u,
        'anom_v': anom_v,
        'clim_u': clim_u,
        'clim_v': clim_v
    }

# Calcular anomalias
results = calculate_september_wind_anomalies(da_u, da_v)

# Selecionar anos específicos em 200 hPa
target_dates = [
    '1982-09-01T00:00:00.000000000',
    '1997-09-01T00:00:00.000000000', 
    '2015-09-01T00:00:00.000000000',
    '2023-09-01T00:00:00.000000000'
]

anoms = [results['anom_u'].sel(valid_time=date, pressure_level=200) for date in target_dates]



fig, axs = plt.subplots(2, 2, figsize=(16, 8), 
                       subplot_kw={'projection': ccrs.Orthographic(-50, -60)})

# Níveis de contorno otimizados para vento em 200, 500 e 850 hPa
levels = np.linspace(-15, 15, 13)
cmap = plt.get_cmap('jet').copy()
cmap.set_over('darkred')
cmap.set_under('darkblue')

for ax, data, date in zip(axs.flat, anoms, target_dates):
    year = pd.to_datetime(date).year
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    
    # Plot com foco na América do Sul
    cf = ax.contourf(da_v.longitude, da_v.latitude, data,
                    levels=levels, cmap=cmap,
                    extend='both', transform=ccrs.PlateCarree())
    
    ax.set_title(f'Setembro {year}', fontsize=16)
    ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('m s-¹', fontsize=20)
cbar.ax.tick_params(labelsize=20)

#plt.suptitle('Anomalias de vento meridional em Setembro - Climatologia 1991-2020 (ERA5)',
 #           fontsize=14, y=1.02)
plt.tight_layout()
plt.show()
