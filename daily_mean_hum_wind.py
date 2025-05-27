import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import pandas as pd
import cdsapi
import zipfile

dataset = "derived-era5-pressure-levels-daily-statistics"
request = {
    "product_type": "reanalysis",
    "variable": [
        "specific_humidity",
        "u_component_of_wind",
        "v_component_of_wind"
    ],
    "year": "2023",
    "month": ["09"],
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
        "28", "29", "30"
    ],
    "pressure_level": ["850"],
    "daily_statistic": "daily_mean",
    "time_zone": "utc+03:00",
    "frequency": "1_hourly"
}
url= 'https://cds.climate.copernicus.eu/api'
key= 'a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a'
client = cdsapi.Client(url,key)
client.retrieve(dataset, request).download()

with zipfile.ZipFile('91daed333dcb1c64ac782b61389720e9.zip', 'r') as zip_ref:
    zip_ref.extractall('derived-era5-pressure-levels-daily-statistics')
    arquivos = zip_ref.namelist()
    print(arquivos)
    
ds_u = ('derived-era5-pressure-levels-daily-statistics/u_component_of_wind_0_daily-mean.nc')
ds_v = ('derived-era5-pressure-levels-daily-statistics/v_component_of_wind_0_daily-mean.nc')
ds_sp =('derived-era5-pressure-levels-daily-statistics/specific_humidity_stream-oper_daily-mean.nc')

ds_open_u = xr.open_dataset(ds_u)
ds_open_v = xr.open_dataset(ds_v)
ds_open_sp = xr.open_dataset(ds_sp)

ds_180_u = ds_open_u.assign_coords(longitude=(((ds_open_u.longitude + 180) % 360) - 180)).sortby('longitude')
ds_180_v = ds_open_v.assign_coords(longitude=(((ds_open_v.longitude + 180) % 360) - 180)).sortby('longitude')
ds_180_sp = ds_open_sp.assign_coords(longitude=(((ds_open_sp.longitude + 180) % 360) - 180)).sortby('longitude')


da_u = ds_180_u['u']
da_v = ds_180_v['v']
da_q = ds_180_sp['q']

def kg_to_g(q):
    return q * 1000

hum_g = kg_to_g(da_q)

wind_magnitude = np.sqrt(da_u**2 + da_v**2)

data = '2023-09-14'
wind_selected = wind_magnitude.sel(valid_time=data, pressure_level=850)
q_selected = hum_g.sel(valid_time=data, pressure_level=850)

magnitude_levels = np.linspace(0, 16, 13)

U_selected = da_u.sel(valid_time=data, pressure_level=850)
V_selected = da_v.sel(valid_time=data, pressure_level=850)

mask = wind_selected > 12
U_filtrado = U_selected.where(mask)
V_filtrado = V_selected.where(mask)

lon = wind_magnitude.longitude
lat = wind_magnitude.latitude

fig, ax = plt.subplots(figsize=(16,8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(cfeature.COASTLINE)
ax.set_extent([-20, -80, -50, 0], crs=ccrs.PlateCarree())

hum_contour = ax.contourf(da_q.longitude, da_q.latitude, q_selected, levels=magnitude_levels,
                           cmap='Blues', extend='max', transform=ccrs.PlateCarree())

pular = 6
ax.quiver(
    lon[::pular], lat[::pular],
    U_filtrado[::pular, ::pular].values,
    V_filtrado[::pular, ::pular].values,
    transform=ccrs.PlateCarree(), scale=600, color='black',
    headwidth=2, headlength=3, linewidth=8, alpha=1
)
ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False),
              linewidth=1, color='gray', alpha=0.5, linestyle='--', transform=ccrs.PlateCarree())

colorbar_ticks = np.linspace(0,16,7)
colorbar = plt.colorbar(hum_contour, ticks=colorbar_ticks, orientation='horizontal', pad=0.1, aspect=40, shrink=0.7)
colorbar.set_label('g kg-¹',size=18)