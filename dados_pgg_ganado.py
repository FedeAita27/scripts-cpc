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
import esmtools.stats as esm
import matplotlib.gridspec as gridspec
from matplotlib import colors



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
    "month": ["01", "02", "03", "04", "05", "06",
              "07", "08", "09", "10", "11", "12"],
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [-20, -60, -36, -45]
}

url='https://cds.climate.copernicus.eu/api'
key='a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a'
client = cdsapi.Client(url,key)
client.retrieve(dataset, request).download()

ds_3 = ('bcae191c3fb7d9779d19b2e7f97d1308.nc')

ds_open_3 = xr.open_dataset(ds_3) # swvl1, swvl2, t2m e tp
ds_open_3

ds_180 = ds_open_3.assign_coords(longitude=(((ds_open_3.longitude + 180) % 360) - 180)).sortby('longitude')
da_var = ds_180['t2m']
da_var
da_var[0,:,:].plot()

def k_c(K):
    return K - 273.15

c = k_c(da_var)
da_var=c
'''
def m_mm(m):
    return m * 1000
tp_mm = m_mm(da_var)
da_var = tp_mm
 '''   
def calcular_anomalia_decadal_anual_detrended_deseasonalized(da_var, climatology_period=(1991, 2020)):
    # Remover sazonalidade: subtrair a média mensal climatológica (média de cada mês)
    da_var.coords['month'] = da_var['valid_time'].dt.month
    media_mensal = da_var.groupby('month').mean('valid_time')
    deseasonalized = da_var.groupby('month') - media_mensal

    # Extrair ano
    deseasonalized.coords['year'] = deseasonalized['valid_time'].dt.year

    # Média anual da série dessazonalizada
    anual = deseasonalized.groupby('year').mean('valid_time')

    # Criar série auxiliar de tempo (1, 2, ..., N)
    tempo = anual.year.values
    var_aux = xr.Dataset({'serie': ('year', np.arange(1, len(tempo)+1)), 'year': tempo})

    # Regressão linear para remover tendência
    reg = esm.linregress(x=var_aux, y=anual, dim='year')

    # Acessar parâmetros da regressão via índice (slope: 0, intercept: 1)
    anual_detrended = anual - (reg.isel(parameter=0) * var_aux['serie'] + reg.isel(parameter=1))

    # Climatologia anual sem tendência e sazonalidade
    climatologia = anual_detrended.sel(year=slice(*climatology_period)).mean('year')

    # Anomalia anual
    anomalia_ano = anual_detrended - climatologia

    # Anomalias decadais
    decadas = {
        '1995-2004': slice(1995, 2004),
        '2005-2014': slice(2005, 2014),
        '2015-2024': slice(2015, 2024)
    }

    anomalias_decadais = {
        dec: anomalia_ano.sel(year=intervalo).mean('year')
        for dec, intervalo in decadas.items()
    }

    return {
        'anomalia_ano': anomalia_ano,
        'anomalias_decadais': anomalias_decadais,
        'climatologia': climatologia,
        'anual_detrended': anual_detrended
    }

resultados = calcular_anomalia_decadal_anual_detrended_deseasonalized(da_var)

# Lista de décadas
decadas = ['1995-2004', '2005-2014', '2015-2024']

# Feição de contorno de estados
reg = cfeature.NaturalEarthFeature(category='cultural',
                                   name='admin_1_states_provinces_lines',
                                   scale='50m', facecolor='none')

# Criar figura e subplots
fig, axs = plt.subplots(1, 3, figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
levels = np.linspace(-2, 2, 11)  # escala de anomalia T2M

for i, decada in enumerate(decadas):
    da_var_selected = resultados['anomalias_decadais'][decada]
    ax = axs[i]

    # Adicionar feições
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
    ax.add_feature(reg)

    ax.set_extent([-60, -45, -35, -20], crs=ccrs.PlateCarree())

    # Plotar anomalia
    cf = ax.contourf(da_var_selected.longitude, da_var_selected.latitude, da_var_selected,
                     levels=levels, cmap='viridis_r',
                     extend='both', transform=ccrs.PlateCarree())

    # Grid
    gl = ax.gridlines(draw_labels=True, linewidth=0.5,
                      linestyle='--', color='gray', alpha=0.5)
    gl.top_labels = False
    gl.bottom_labels = True
    gl.left_labels = (i == 0)
    gl.right_labels = False

    # Título
    ax.set_title(f'Anomalia T2M\n{decada}', fontsize=16, color='black')

# Barra de cores comum
cbar = fig.colorbar(cf, ax=axs, orientation='horizontal', pad=0.05, shrink=0.8)
cbar.set_label('Anomalia T2M (°C)', fontsize=20)
cbar.ax.tick_params(labelsize=15)

plt.tight_layout()
plt.show()

'''
fig, ax = plt.subplots(figsize=(16,8), subplot_kw={'projection': ccrs.PlateCarree()})

ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
ax.add_feature(reg)
ax.set_extent([-60, -45, -35, -20], crs=ccrs.PlateCarree())
levels= np.linspace(0,0.5,6)

cf = ax.contourf(da_var_selected.longitude, da_var_selected.latitude, da_var_selected,
                    levels=levels, cmap='viridis_r',
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
'''
plt.show()
