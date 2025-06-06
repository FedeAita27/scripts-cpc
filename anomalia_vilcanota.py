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
import cdsapi
from scipy.stats import linregress
import statsmodels.api as sm
from statsmodels.tsa.seasonal import STL
import rasterio
from rasterio.plot import show
import geopandas as gpd

dataset = "reanalysis-era5-single-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": ["2m_temperature"],
    "year": [
        "1991", "1992", "1993",
        "1994", "1995", "1996",
        "1997", "1998", "1999",
        "2000", "2001", "2002",
        "2003", "2004", "2005",
        "2006", "2007", "2008",
        "2009", "2010", "2011",
        "2012", "2013", "2014",
        "2015", "2016", "2017",
        "2018", "2019", "2020",
        "2021", "2022", "2023",
        "2024", "2025"
    ],
    "month": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12"
    ],
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [-13, -72, -15, -70]
}

url='https://cds.climate.copernicus.eu/api'
key='a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a'
client = cdsapi.Client(url,key)
client.retrieve(dataset, request).download()

nc = 'b2137d03f35aaa3ed2ee04062ab4bebb.nc'

ds_t925 = xr.open_dataset(nc)
 
ds_180 = ds_t925.assign_coords(longitude=(((ds_t925.longitude + 180) % 360) - 180)).sortby('longitude')
da_t925 = ds_180['t2m']

da_degc =  da_t925 - 273.15
da_t = da_degc

da_t[0,:,:].plot()

def calculate_t2m_anomalies(da_t, reference_period=('1991-01-01', '2020-12-01')):
    # 1. Filtrar apenas dados de setembro
    sept_t = da_t.sel(valid_time=da_t['valid_time.month'] == 9)
    
    # 2. Calcular climatologia (média 1991-2020)
    clim_t = sept_t.sel(valid_time=slice(reference_period[0], reference_period[1])).groupby('valid_time.month').mean()
    
    # 3. Calcular anomalias (dados observados - climatologia)
    anom_t = sept_t.groupby('valid_time.month') - clim_t
    
    return {
        'anom_t': anom_t,
        'clim_t': clim_t,
    }

# Calcular anomalias
results = calculate_t2m_anomalies(da_t)

# Selecionar anos específicos em 200 hPa
#target_date = '2023-09-01T00:00:00.000000000'

da_t_selected = da_t.sel(valid_time='2024-01-01')

srtm_path = 'C:/Users/Aluno/Desktop/Federico/BIC_IC -  Kátia/s14_w071_1arc_v3.tif'
with rasterio.open(srtm_path) as src:
    fig, ax = plt.subplots(figsize=(10, 8))
    show(src, ax=ax, cmap='terrain')


    # Opcional: sobrepor com shapefile (como limites administrativos)
    shapefile_path =  'C:/Users/Aluno/Desktop/Federico/BIC_IC -  Kátia/bacia - CopiaUTM.shp'

    gdf = gpd.read_file(shapefile_path)
    gdf.boundary.plot(ax=ax, color='black', linewidth=0.5)

    ax.set_title('Modelo Digital de Elevação - SRTM')
    plt.show()


#anom = results['anom_t925'].sel(valid_time=target_date, pressure_level=925)

fig, ax = plt.subplots(figsize=(16, 8), 
                       subplot_kw={'projection': ccrs.PlateCarree()})

# Níveis de contorno otimizados para vento em 200hPa
levels = np.linspace(0, 25, 13)
cmap = plt.get_cmap('coolwarm').copy()

ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    
    # Plot com foco na América do Sul
cf = ax.contourf(da_t.longitude, da_t.latitude, da_t_selected,
                    levels=levels, cmap=cmap,
                    extend='both', transform=ccrs.PlateCarree())

#year = pd.to_datetime(target_date).year    
#ax.set_title(f'Setembro {year}', fontsize=16)
ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5)

# Barra de cores
#cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, orientation='vertical', extend='both')
cbar.set_label('°C', fontsize=20)
cbar.ax.tick_params(labelsize=20)

#plt.suptitle('Anomalias de vento meridional em Setembro - Climatologia 1991-2020 (ERA5)',
 #           fontsize=14, y=1.02)
plt.tight_layout()
plt.show()
