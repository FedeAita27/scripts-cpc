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

########## ANOMALIA TCWV #########
tcwvnetcdf4 = ('C:/Users/feder/OneDrive/Documentos/nc/sphum_september.nc')

ds_tcwv = xr.open_dataset(tcwvnetcdf4)
ds_tcwv

ds_180_tcwv = ds_tcwv.assign_coords(longitude=(((ds_tcwv.longitude + 180) % 360) - 180)).sortby('longitude')
da_tcwv = ds_180_tcwv['q']
da_tcwv

def g_kg(q):
    q / 1000.0
    return q
g_kg = g_kg(da_tcwv)
da_tcwv = g_kg

da_tcwv[0,:,:].plot()

def calculate_september_q_anomalies(da_tcwv, reference_period=('1991-09-01', '2020-09-01')):
    # 1. Filtrar apenas dados de setembro
    sept_tcwv = da_tcwv.sel(valid_time=da_tcwv['valid_time.month'] == 9)
    
    # 2. Calcular climatologia (média 1991-2020)
    clim_tcwv = sept_tcwv.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    
    # 3. Calcular anomalias (dados observados - climatologia)
    anom_tcwv = sept_tcwv.groupby('valid_time.month') - clim_tcwv
    
    return {
        'anom_tcwv': anom_tcwv,
        'clim_tcwv': clim_tcwv,
    }

######### ANOMALIA VENTO #########
windsnetcdf4 = ('C:/Users/feder/OneDrive/Documentos/nc/winds_september.nc')

ds_wind = xr.open_dataset(windsnetcdf4)
 
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
 
    anom_wind = np.sqrt(anom_u**2 + anom_v**2)
    
    return {
        'anom_u': anom_u,
        'anom_v': anom_v,
        'clim_u': clim_u,
        'clim_v': clim_v,
        'anom_wind': anom_wind
    }

# Calcular anomalias
results_winds = calculate_september_wind_anomalies(da_u, da_v)

# Selecionar anos específicos em 200 hPa
target_date = '2023-09-01T00:00:00.000000000'

anom_w = results_winds['anom_wind'].sel(valid_time=target_date, pressure_level=850)

# Calcular anomalias
results_tcwv = calculate_september_q_anomalies(da_tcwv)

# Selecionar anos específicos em 200 hPa
anom_hum = results_tcwv['anom_tcwv'].sel(valid_time=target_date)

mask = anom_w > 4
U_filtered = anom_w.where(mask)
V_filtered = anom_w.where(mask)

fig, ax = plt.subplots(figsize=(16, 8), 
                       subplot_kw={'projection': ccrs.Orthographic(-50, -60)})

# Níveis de contorno otimizados
levels = np.linspace(-9,9, 13)
cmap = plt.get_cmap('Blues').copy()


ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    
    # Plot com foco na América do Sul
cf = ax.contourf(da_tcwv.longitude, da_tcwv.latitude, anom_hum,
                    levels=levels, cmap=cmap,
                    extend='both', transform=ccrs.PlateCarree())
pular = 8
vet = ax.quiver(ds_wind.longitude[::10], ds_wind.latitude[::10],
          U_filtered.values[::pular, ::pular], V_filtered.values[::pular, ::pular],
          color='black', scale=600, headwidth=2, headlength=3, linewidth=1, alpha=1,
          transform=ccrs.PlateCarree())

year = pd.to_datetime(target_date).year    
#ax.set_title(f'Setembro {year}', fontsize=16)
ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal', extend='both')
cbar.set_label('mm', fontsize=20)
cbar.ax.tick_params(labelsize=20)

#plt.suptitle('Anomalias de vento meridional em Setembro - Climatologia 1991-2020 (ERA5)',
 #           fontsize=14, y=1.02)
plt.tight_layout()
plt.show()

