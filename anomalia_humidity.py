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

humnetcdf4 = ('C:/Users/feder/OneDrive/Documentos/nc/sphum_september.nc')

ds_hum = xr.open_dataset(humnetcdf4)

ds_180 = ds_hum.assign_coords(longitude=(((ds_hum.longitude + 180) % 360) - 180)).sortby('longitude')
da_hum = ds_180['q']

def kg_to_g(q):
    return q * 1.0

hum_g = kg_to_g(da_hum)


hum_g[0,:,:].plot()

def calculate_september_q_anomalies(hum_g, reference_period=('1991-09-01', '2020-09-01')):
    # 1. Filtrar apenas dados de setembro
    sept_q = hum_g.sel(valid_time=hum_g['valid_time.month'] == 9)
    
    # 2. Calcular climatologia (média 1991-2020)
    clim_q = sept_q.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    
    # 3. Calcular anomalias (dados observados - climatologia)
    anom_q = sept_q.groupby('valid_time.month') - clim_q
    
    return {
        'anom_q': anom_q,
        'clim_q': clim_q,
    }

# Calcular anomalias
results = calculate_september_q_anomalies(hum_g)

# Selecionar anos específicos em 200 hPa
target_dates = [
    '1982-09-01T00:00:00.000000000',
    '1997-09-01T00:00:00.000000000', 
    '2015-09-01T00:00:00.000000000',
    '2023-09-01T00:00:00.000000000'
]

anoms = [results['anom_q'].sel(valid_time=date, pressure_level=850) for date in target_dates]

minim = min(hum_g)
maxim = max(hum_g)

fig, axs = plt.subplots(2, 2, figsize=(16, 8), 
                       subplot_kw={'projection': ccrs.Orthographic(-50, -60)})

# Níveis de contorno otimizados
levels = np.linspace(0,15, 13)
cmap = plt.get_cmap('Blues').copy()
cmap.set_over('darkred')
cmap.set_under('darkblue')

for ax, data, date in zip(axs.flat, anoms, target_dates):
    year = pd.to_datetime(date).year
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    
    # Plot com foco na América do Sul
    cf = ax.contourf(hum_g.longitude, hum_g.latitude, data,
                    levels=levels, cmap=cmap,
                    extend='both', transform=ccrs.PlateCarree())
    
    ax.set_title(f'Setembro {year}', fontsize=12)
    ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('g kg^-1', fontsize=12)

plt.suptitle('Anomalia de umidade específica em setembro - Climatologia 1991-2020 (ERA5)',
            fontsize=14, y=1.02)
plt.tight_layout()
plt.show()
