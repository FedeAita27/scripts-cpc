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

mslnetcdf4 = ('C:/Users/feder/OneDrive/ERA5/TCC/nc/msl_september.nc')

ds_msl = xr.open_dataset(mslnetcdf4)
ds_msl
 
ds_180 = ds_msl.assign_coords(longitude=(((ds_msl.longitude + 180) % 360) - 180)).sortby('longitude')
da_msl = ds_180['msl']

def Pa_to_mb(Pa):
   return Pa / 100.0

msl_mb = Pa_to_mb(da_msl)
da_msl = msl_mb

da_msl[0,:,:].plot()

def calculate_september_msl_anomalies(da_msl, reference_period=('1991-09-01', '2020-09-01')):
    # 1. Filtrar apenas dados de setembro
    sept_msl = da_msl.sel(valid_time=da_msl['valid_time.month'] == 9)
    
    # 2. Calcular climatologia (média 1991-2020)
    clim_msl = sept_msl.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    
    # 3. Calcular anomalias (dados observados - climatologia)
    anom_msl = sept_msl.groupby('valid_time.month') - clim_msl
    
    return {
        'anom_msl': anom_msl,
        'clim_msl': clim_msl,
    }

# Calcular anomalias
results = calculate_september_msl_anomalies(da_msl)

# Selecionar anos específicos em 200 hPa
target_dates = [
    '1982-09-01T00:00:00.000000000',
    '1997-09-01T00:00:00.000000000', 
    '2015-09-01T00:00:00.000000000',
    '2023-09-01T00:00:00.000000000'
]

anoms = [results['anom_msl'].sel(valid_time=date) for date in target_dates]

fig, axs = plt.subplots(2, 2, figsize=(16, 8), 
                       subplot_kw={'projection': ccrs.Orthographic(-50, -60)})

# Níveis de contorno otimizados para vento em 200hPa
levels = np.linspace(-10, 10, 13)
cmap = plt.get_cmap('coolwarm').copy()
#cmap.set_over('darkred')
#cmap.set_under('darkblue')

for ax, data, date in zip(axs.flat, anoms, target_dates):
    year = pd.to_datetime(date).year
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    
    # Plot com foco na América do Sul
    cf = ax.contourf(da_msl.longitude, da_msl.latitude, data,
                    levels=levels, cmap=cmap,
                    extend='both', transform=ccrs.PlateCarree())
    
    ax.set_title(f'Setembro {year}', fontsize=16, color='black')
    ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('hPa', fontsize=20)
cbar.ax.tick_params(labelsize=20)

#plt.suptitle('Anomalias pressão atmosférica ao nível do mar em Setembro - Climatologia 1991-2020 (ERA5)',
 #           fontsize=14, y=1.02)
plt.tight_layout()
plt.show()
