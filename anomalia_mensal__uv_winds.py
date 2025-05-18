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

windsnetcdf4 = ('C:/Users/Aluno/Downloads/6b2530fe14aef886a7e97a22d2f1dc93.nc')

ds_wind = xr.open_dataset(windsnetcdf4)
ds_wind
#ds_wind['date']
 
ds_180 = ds_wind.assign_coords(longitude=(((ds_wind.longitude + 180) % 360) - 180)).sortby('longitude')
da_v = ds_180['v']
da_u = ds_180['u']

da_v[0,:,:].plot()
def calculate_wind_anomalies(da_u, da_v, reference_period=('1991-01-01', '2020-12-31')):
    da_v['valid_time'] = pd.to_datetime(da_v['valid_time'].values, format='%Y%m%d')
    da_u['valid_time'] = pd.to_datetime(da_u['valid_time'].values, format='%Y%m%d')

    yearly_mean_V = da_v.groupby('valid_time.year').mean(keep_attrs=True)
    yearly_mean_U = da_u.groupby('valid_time.year').mean(keep_attrs=True)
    
    ref_winds_V = yearly_mean_V.where((yearly_mean_V.year > 1990) & (yearly_mean_V.year < 2021), drop=True)
    ref_winds_U = yearly_mean_U.where((yearly_mean_U.year > 1990) & (yearly_mean_U.year < 2021), drop=True)

    ref_mean_winds_V = ref_winds_V.mean(dim="year", keep_attrs=True)
    ref_mean_winds_U = ref_winds_U.mean(dim="year", keep_attrs=True)

    clim_period_V = da_v.sel(valid_time=slice('1991-01-01', '2020-12-31'))
    clim_period_U = da_u.sel(valid_time=slice('1991-01-01', '2020-12-31'))

    clim_month_V = clim_period_V.groupby('valid_time.month').mean()
    clim_month_U = clim_period_U.groupby('valid_time.month').mean()

    weights_V = np.cos(np.deg2rad(clim_month_V.latitude))
    weights_U = np.cos(np.deg2rad(clim_month_U.latitude))

    weights_V.name = "weights_V"
    weights_U.name = "weights_U"

    clim_month_weighted_V = clim_month_V.weighted(weights_V)
    clim_month_weighted_U = clim_month_U.weighted(weights_U)
    mean_V = clim_month_weighted_V.mean(["longitude", "latitude"])
    mean_U = clim_month_weighted_U.mean(["longitude", "latitude"])

    anom_month_V = da_v.groupby('valid_time.month') - clim_month_V
    anom_month_U = da_u.groupby('valid_time.month') - clim_month_U
    
    return {
        'anom_month_U': anom_month_U,
        'anom_month_V': anom_month_V,
        'clim_month_U': clim_month_U,
        'clim_month_V': clim_month_V
    }

results = calculate_wind_anomalies(da_v, da_u)

anom_month_U = results["anom_month_U"]
anom_month_V = results["anom_month_V"]
clim_month_U = results["clim_month_U"]
clim_month_V = results["clim_month_V"]

anom_1 = anom_month_U.sel(valid_time='2016-09-01', pressure_level=200)
#anom_2 = anom_month_U.sel(date='2023-11-01', pressure_level=200)
#anom_3 = anom_month_U.sel(date='2024-04-01', pressure_level=200)

fig, axs = plt.subplots(1, 1, figsize = (16, 8), subplot_kw={'projection': ccrs.Orthographic(-50,-60)})

def plot_anomaly(ax, data, lon, lat, title):
    magnitude_levels = np.linspace(-18, 18, 8)
    ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
    ax.add_feature(cfeature.COASTLINE, linewidth=1)
    ax.set_title(title, size=10)
    wind_contour = ax.contourf(lon, lat, data, cmap='coolwarm',
                 levels=magnitude_levels,
                 extend='both', transform=ccrs.PlateCarree())

    ax.gridlines(draw_labels=dict(left=False, bottom=False, top=False, right=False),
              linewidth=1, color='gray', alpha=0.5, linestyle='--', transform=ccrs.Orthographic())
    colorbar_ticks = np.linspace(-18, 18, 8)
    cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.05])  # [left, bottom, width, height]
    cbar = fig.colorbar(wind_contour, cax=cbar_ax, orientation='horizontal', ticks=colorbar_ticks)
    cbar.set_label('Magnitude da Anomalia de Vento (m/s)', fontsize=12)

plot_level = 200

# OLHAR AS LINHAS 71, 72 e 73
reference_month_1 = 'setembro'
#reference_month_2 = 'novembro'
#reference_month_3 = 'abril'

plot_anomaly(axs[0], anom_1, da_u.longitude, da_u.latitude, f'Anomalia de vetor vento em {plot_level} hPa - {reference_month_1} 2023')
#plot_anomaly(axs[1], anom_2, da_u.longitude, da_u.latitude, f'Anomalia de vetor vento em {plot_level} hPa - {reference_month_2} 2023')
#plot_anomaly(axs[2], anom_3, da_u.longitude, da_u.latitude, f'Anomalia de vetor vento em {plot_level} hPa - {reference_month_3} 2024')

#colorbar_ticks = np.linspace(-6, 6, 10)
#cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.1])  # [left, bottom, width, height]
#cbar = fig.colorbar(wind_contour, cax=cbar_ax, orientation='horizontal', ticks=colorbar_ticks)
#cbar.set_label('Magnitude da Anomalia de Vento (m/s)')

#plt.subplots_adjust(wspace=0.3, bottom=0.5)
