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

tsmnetcdf4 = ('C:/Users/feder/OneDrive/ERA5/TCC/nc/tsm_september.nc')

ds_tsm = xr.open_dataset(tsmnetcdf4)
ds_tsm
 
ds_180 = ds_tsm.assign_coords(longitude=(((ds_tsm.longitude + 180) % 360) - 180)).sortby('longitude')
da_tsm = ds_180['sst']

da_degc =  da_tsm - 273.15
da_tsm = da_degc

da_tsm[0,:,:].plot()
def calculate_september_tsm_anomalies(da_tsm, reference_period=('1991-09-01', '2020-09-01')):
    # 1. Filtrar apenas dados de setembro
    sept_tsm = da_tsm.sel(valid_time=da_tsm['valid_time.month'] == 9)
    
    # 2. Calcular climatologia (média 1991-2020)
    clim_tsm = sept_tsm.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    
    # 3. Calcular anomalias (dados observados - climatologia)
    anom_tsm = sept_tsm.groupby('valid_time.month') - clim_tsm
    
    return {
        'anom_tsm': anom_tsm,
        'clim_tsm': clim_tsm
    }

# Calcular anomalias
results = calculate_september_tsm_anomalies(da_tsm)

# Selecionar anos específicos em 200 hPa
target_dates = [
    #'1997-09-01T00:00:00.000000000', 
    #'2015-09-01T00:00:00.000000000',
    '2023-09-01T00:00:00.000000000'
]

anoms = [results['anom_tsm'].sel(valid_time=date) for date in target_dates]

fig, ax = plt.subplots(1, 1, figsize=(16, 8), 
                   subplot_kw={'projection': ccrs.PlateCarree()})

# Níveis de contorno otimizados para vento em 200hPa
levels = np.linspace(-6, 6, 13)
cmap = plt.get_cmap('coolwarm').copy()
#cmap.set_over('darkred')
#cmap.set_under('darkblue')

# Como agora temos apenas um eixo, não precisamos do loop
ax.stock_img()
year = pd.to_datetime(target_dates[0]).year  # assumindo que target_dates tem pelo menos um elemento
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)

# Plot com foco na América do Sul
cf = ax.contourf(da_tsm.longitude, da_tsm.latitude, anoms[0],  # assumindo que anoms tem pelo menos um elemento
                levels=levels, cmap=cmap,
                extend='both', transform=ccrs.PlateCarree())

#ax.set_title(f'Setembro {year}', fontsize=12)
ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False), linewidth=0.5, color='gray', alpha=0.5)
cbar_width = 1
cbar_left = 0.5 - (cbar_width / 2)

# Barra de cores
cbar_ax = fig.add_axes([cbar_left, -0.05, cbar_width, 0.04])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('°C', fontsize=22)
cbar.ax.tick_params(labelsize=22)

#plt.suptitle('Anomalia de TSM em Setembro de 2023 - Climatologia 1991-2020 (ERA5)',
 #           fontsize=20)
plt.tight_layout()
plt.show()

'''
fig, axs = plt.subplots(1, 1, figsize=(16, 8), 
                       subplot_kw={'projection': ccrs.PlateCarree()})

# Níveis de contorno otimizados para vento em 200hPa
levels = np.linspace(-6, 6, 13)
cmap = plt.get_cmap('coolwarm').copy()
cmap.set_over('darkred')
cmap.set_under('darkblue')

for ax, data, date in zip(axs.flat, anoms, target_dates):
    ax.stock_img()
    year = pd.to_datetime(date).year
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    
    # Plot com foco na América do Sul
    cf = ax.contourf(da_tsm.longitude, da_tsm.latitude, data,
                    levels=levels, cmap=cmap,
                    extend='both', transform=ccrs.PlateCarree())
    
    ax.set_title(f'Setembro {year}', fontsize=12)
    ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('°C', fontsize=12)

plt.suptitle('Anomalia TSM em Setembro - Climatologia 1991-2020 (ERA5)',
            fontsize=20, y=1.02)
plt.tight_layout()
plt.show()
'''