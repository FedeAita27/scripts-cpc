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

gptnetcdf4 = ('C:/Users/feder/OneDrive/Documentos/nc/gpt_setembro.nc')

ds_gpt = xr.open_dataset(gptnetcdf4)
ds_gpt
 
ds_180 = ds_gpt.assign_coords(longitude=(((ds_gpt.longitude + 180) % 360) - 180)).sortby('longitude')
da_gpt = ds_180['z']

def z_g(z):
   return z / 9.80665

z_g = z_g(da_gpt)
da_gpt = z_g

da_gpt[0,:,:].plot()

def calculate_september_gpt_anomalies(da_gpt, reference_period=('1991-09-01', '2020-09-01')):
    # 1. Filtrar apenas dados de setembro
    sept_gpt = da_gpt.sel(valid_time=da_gpt['valid_time.month'] == 9)
    
    # 2. Calcular climatologia (média 1991-2020)
    clim_gpt = sept_gpt.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    
    # 3. Calcular anomalias (dados observados - climatologia)
    anom_gpt = sept_gpt.groupby('valid_time.month') - clim_gpt
    
    return {
        'anom_gpt': anom_gpt,
        'clim_gpt': clim_gpt,
    }

# Calcular anomalias
results = calculate_september_gpt_anomalies(da_gpt)

# Selecionar anos específicos em 200 hPa
target_date = '2023-09-01T00:00:00.000000000'

anom = results['anom_gpt'].sel(valid_time=target_date, pressure_level=200)

fig, ax = plt.subplots(figsize=(16, 8), 
                       subplot_kw={'projection': ccrs.Orthographic(-50, -60)})

# Níveis de contorno otimizados para vento em 200hPa
levels = np.linspace(-120, 120, 13)
cmap = plt.get_cmap('coolwarm').copy()

year = pd.to_datetime(target_date).year
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    
# Plot com foco na América do Sul
cf = ax.contourf(da_gpt.longitude, da_gpt.latitude, anom,
                    levels=levels, cmap=cmap,
                    extend='both', transform=ccrs.PlateCarree())
    
ax.set_title(f'Setembro {year}', fontsize=16, color='black')
ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('mm', fontsize=20)
cbar.ax.tick_params(labelsize=20)

#plt.suptitle('Anomalias pressão atmosférica ao nível do mar em Setembro - Climatologia 1991-2020 (ERA5)',
 #           fontsize=14, y=1.02)
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(16, 8), 
                       subplot_kw={'projection': ccrs.PlateCarree(-50, -30)})

# Níveis de contorno otimizados para vento em 200hPa
levels = np.linspace(-120, 120, 13)
cmap = plt.get_cmap('coolwarm').copy()

year = pd.to_datetime(target_date).year
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    
# Plot com foco na América do Sul
cf = ax.contourf(da_gpt.longitude, da_gpt.latitude, anom,
                    levels=levels, cmap=cmap,
                    extend='both', transform=ccrs.PlateCarree())
    
ax.set_title(f'Setembro {year}', fontsize=16, color='black')
ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('mm', fontsize=20)
cbar.ax.tick_params(labelsize=20)

#plt.suptitle('Anomalias pressão atmosférica ao nível do mar em Setembro - Climatologia 1991-2020 (ERA5)',
 #           fontsize=14, y=1.02)
plt.tight_layout()
plt.show()