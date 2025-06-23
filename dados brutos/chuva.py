import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import xarray as xr
import numpy as np
import pandas as pd
import scipy
from netCDF4 import Dataset
import cdsapi

dataset = "derived-era5-single-levels-daily-statistics"
request = {
    "product_type": "reanalysis",
    "variable": ["convective_rain_rate"],
    "year": "2024",
    "month": ["04", "05"],
    "day": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12",
        "13", "14", "15",
        "16", "17", "18",
        "19", "20", "21",
        "22", "23", "24",
        "25", "26", "27",
        "28", "29", "30",
        "31"
    ],
    "daily_statistic": "daily_mean",
    "time_zone": "utc+03:00",
    "frequency": "1_hourly"
}

url='https://cds.climate.copernicus.eu/api'
key='a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a'

client = cdsapi.Client(url, key)
client.retrieve(dataset, request).download()


ds = '28b8df6f3d389275e32d30994993a450.nc'
ds_cr = xr.open_dataset(ds)
ds_180 = ds_cr.assign_coords(longitude=(((ds_cr.longitude + 180) % 360) - 180)).sortby('longitude')
da_crr = ds_180['crr']
da_crr[0,:,:].plot()
da_crr = da_crr * 86400
da_crr[0,:,:].plot()



time_desired = '2024-05-01'
crr = da_crr.sel(valid_time=time_desired)
levels = np.linspace(0, 400, 13)

crr_periodo = da_crr.sel(valid_time=slice('2024-04-26', '2024-05-31'))
crr_soma = crr_periodo.sum(dim='valid_time')


levels = np.linspace(0, 500, 11)

# Define a projeção do mapa
fig, ax = plt.subplots(figsize=(11, 6), subplot_kw={'projection': ccrs.PlateCarree()})

img = ax.contourf(crr_soma.longitude, crr_soma.latitude,
                  crr_soma, cmap='YlGnBu', extend = 'neither', 
                  levels=levels, transform=ccrs.PlateCarree()
                  )

ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False))


ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)

ax.set_title('Precipitação maio 2024 - 26/04 a 05/05', fontsize=13)
ax.set_extent([-45, -60, -20, -35], crs=ccrs.PlateCarree())

colorbar_ticks = np.linspace(0, 500, 6)

cbar = plt.colorbar(img, ax=ax, orientation='vertical',
                    ticks=colorbar_ticks, pad=0.05, aspect=20)
cbar.set_label('mm', size=16)

plt.tight_layout()
plt.show()