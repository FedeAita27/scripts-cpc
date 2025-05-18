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


t925netcdf4 = ('C:/Users/feder/OneDrive/Documentos/nc/t925_september.nc')

ds_t925 = xr.open_dataset(t925netcdf4)
 
ds_180 = ds_t925.assign_coords(longitude=(((ds_t925.longitude + 180) % 360) - 180)).sortby('longitude')
da_t925 = ds_180['t']

da_degc =  da_t925 - 273.15
da_t925 = da_degc

def calculate_september_t925_anomalies(da_t925, reference_period=('1991-09-01', '2020-09-01')):
    # 1. Filtrar apenas dados de setembro
    sept_t925 = da_t925.sel(valid_time=da_t925['valid_time.month'] == 9)
    
    # 2. Calcular climatologia (média 1991-2020)
    clim_t925 = sept_t925.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    
    # 3. Calcular anomalias (dados observados - climatologia)
    anom_t925 = sept_t925.groupby('valid_time.month') - clim_t925
    
    return {
        'anom_t925': anom_t925,
        'clim_t925': clim_t925,
    }

# Calcular anomalias
results = calculate_september_t925_anomalies(da_t925)

# Selecionar anos específicos em 200 hPa
target_date = '2023-09-01T00:00:00.000000000'


anom = results['anom_t925'].sel(valid_time=target_date, pressure_level=925)

fig, ax = plt.subplots(figsize=(16, 8), 
                       subplot_kw={'projection': ccrs.Orthographic(-50, -60)})

# Níveis de contorno otimizados para vento em 200hPa
levels = np.linspace(-6, 6, 13)
cmap = plt.get_cmap('coolwarm').copy()

ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    
    # Plot com foco na América do Sul
cf = ax.contourf(da_t925.longitude, da_t925.latitude, anom,
                    levels=levels, cmap=cmap,
                    extend='both', transform=ccrs.PlateCarree())

year = pd.to_datetime(target_date).year    
#ax.set_title(f'Setembro {year}', fontsize=16)
ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('°C', fontsize=20)
cbar.ax.tick_params(labelsize=20)

#plt.suptitle('Anomalias de vento meridional em Setembro - Climatologia 1991-2020 (ERA5)',
 #           fontsize=14, y=1.02)
plt.tight_layout()
plt.show()
