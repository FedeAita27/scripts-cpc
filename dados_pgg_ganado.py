import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import cdsapi
import zipfile

dataset = "reanalysis-era5-land-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "2m_temperature",
        "total_precipitation",
        "volumetric_soil_water_layer_1",
        "volumetric_soil_water_layer_2"
    ],
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
    "month": ["03", "04", "05"],
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [-20, -60, -36, -45]
}

url='https://cds.climate.copernicus.eu/api'
key='a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a'
client = cdsapi.Client(url,key)
client.retrieve(dataset, request).download()

with zipfile.ZipFile("a7e3d3358d073bfa724265a71f19c37b.zip", 'r') as zip_ref:
    zip_ref.extractall("reanalysis-era5-single-levels-monthly-means")
    arquivos = zip_ref.namelist()
    print(arquivos)

ds_0 = ('C:/Users/feder/OneDrive/Documentos/nc/data_0.nc')
ds_1 = ('C:/Users/feder/OneDrive/Documentos/nc/data_1.nc')
ds_2 = ('C:/Users/feder/OneDrive/Documentos/nc/data_2.nc')
ds_3 = ('C:/Users/feder/OneDrive/Documentos/nc/c9af3fdc2aa8598e11f5f68d6b58dfcc.nc')

ds_open_0 = xr.open_dataset(ds_0) # swvl1 e swvl2
ds_open_1 = xr.open_dataset(ds_1) # t2m
ds_open_2 = xr.open_dataset(ds_2) # tp
ds_open_3 = xr.open_dataset(ds_3) # swvl1, swvl2, t2m e tp
ds_open_3

ds_180 = ds_open_3.assign_coords(longitude=(((ds_open_3.longitude + 180) % 360) - 180)).sortby('longitude')
da_var = ds_180['swvl2']
da_var
da_var[0,:,:].plot()
'''
def m_mm(m):
    return m * 1000
tp_mm = m_mm(da_var)
da_var = tp_mm
 '''   

def calculate_seasonal_anomalies(da_var, season_months=[3, 4, 5], reference_period=('1991-03-01', '2020-05-31')):
    # Filtrar os meses da estação (ex: MAM)
    dados = da_var.sel(valid_time=da_var['valid_time'].dt.month.isin(season_months))

    # Criar uma coordenada auxiliar com o ano
    dados.coords['year'] = dados['valid_time'].dt.year

    # Agrupar por ano e tirar a média do trimestre
    seasonal_mean = dados.groupby('year').mean('valid_time')

    # Calcular climatologia (média de 1991 a 2020)
    climatologia = seasonal_mean.sel(year=slice(1991, 2020)).mean('year')

    # Calcular anomalia: cada ano menos a climatologia
    anomalia = seasonal_mean - climatologia

    return {
        'anomalia': anomalia,
        'climatologia': climatologia,
        'seasonal_mean': seasonal_mean
    }


results = calculate_seasonal_anomalies(da_var)
results

# Outono
mam_anom = results['anomalia'].sel(year=2024)
mam_anom

#son = results['anomalia'].sel(valid_time=results['anomalia']['valid_time'].dt.year == 2024) 

reg = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines',
                                   scale='50m', facecolor='none')

fig, ax = plt.subplots(figsize=(16,8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(reg)
ax.set_extent([-60, -45, -35, -20], crs=ccrs.PlateCarree())
levels= np.linspace(-0.12,0.12,6)

cf = ax.contourf(mam_anom.longitude, mam_anom.latitude, mam_anom,
                    levels=levels, cmap='YlGnBu',
                    extend='neither', transform=ccrs.PlateCarree())

    
#ax.set_title('Outono 2016', fontsize=16, color='black')
gl = ax.gridlines(draw_labels=True, linewidth=0.5,
                  linestyle='--', color='gray', alpha=0.5)

gl.top_labels = False
gl.bottom_labels = True
gl.left_labels = True
gl.right_labels = False

# Barra de cores
#cbar_ax = fig.add_axes([0.2, -0.01, 0.6, 0.02])
cbar = fig.colorbar(cf, extend='both')
cbar.set_label('m³/m³ ', fontsize=20)
cbar.ax.tick_params(labelsize=15)

plt.show()
