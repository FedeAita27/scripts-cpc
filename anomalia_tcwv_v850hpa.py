import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import pandas as pd

# Load TCWV data
tcwv_file = 'C:/Users/feder/OneDrive/Documentos/nc/tcwv_september.nc'
ds_tcwv = xr.open_dataset(tcwv_file)
ds_180 = ds_tcwv.assign_coords(longitude=(((ds_tcwv.longitude + 180) % 360) - 180)).sortby('longitude')
da_tcwv = ds_180['tcwv']

# Load wind data
winds_file = 'C:/Users/feder/OneDrive/Documentos/nc/winds_september.nc'
ds_wind = xr.open_dataset(winds_file)
ds_180_wind = ds_wind.assign_coords(longitude=(((ds_wind.longitude + 180) % 360) - 180).sortby('longitude'))
da_u = ds_180_wind['u']
da_v = ds_180_wind['v']

# Function to calculate September anomalies
def calculate_september_anomalies(da, reference_period=('1991-09-01', '2020-09-01')):
    sept_data = da.sel(valid_time=da['valid_time.month'] == 9)
    clim_data = sept_data.sel(valid_time=slice(*reference_period)).groupby('valid_time.month').mean()
    anom_data = sept_data.groupby('valid_time.month') - clim_data
    return anom_data

# Calculate anomalies
tcwv_anom = calculate_september_anomalies(da_tcwv)
u_anom = calculate_september_anomalies(da_u.sel(pressure_level=850))
v_anom = calculate_september_anomalies(da_v.sel(pressure_level=850))
wind_anom = np.sqrt(u_anom**2+v_anom**2)

# Select target date
target_date = '2023-09-01T00:00:00.000000000'
tcwv_plot = tcwv_anom.sel(valid_time=target_date)
wind_plot = wind_anom.sel(valid_time=target_date)

# Create mask where TCWV anomaly is significant (absolute value > 1mm)
tcwv_mask = np.abs(tcwv_plot) > 1

# Create figure
fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(-50, -60))

# Add map features
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.3)

# Plot TCWV anomalies
levels = np.linspace(-9, 9, 13)
cmap = plt.get_cmap('RdBu')
cf = ax.contourf(tcwv_plot.longitude, tcwv_plot.latitude, tcwv_plot,
                 levels=levels, cmap=cmap, extend='both', 
                 transform=ccrs.PlateCarree())

# Plot wind vectors where TCWV anomaly is significant
quiver_skip = 10  # Skip every N points for clearer visualization
ax.quiver(wind_anom.longitude.values[::quiver_skip], 
          wind_anom.latitude.values[::quiver_skip],
          wind_plot.values[::quiver_skip, ::quiver_skip],
          color='black', scale=400, headwidth=2, headlength=3, linewidth=1, alpha=1, 
          transform=ccrs.PlateCarree())

# Add colorbar
cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('mm', fontsize=20)
cbar.ax.tick_params(labelsize=20)

# Add title
year = pd.to_datetime(target_date).year
#plt.title(f'September {year} - TCWV Anomalies and 850 hPa Wind Anomalies', 
  #        fontsize=14, pad=20)

plt.tight_layout()
plt.show()