!pip install geopandas

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
import geopandas as gpd


tpnetcdf4 = Dataset('C:/Users/feder/Downloads/1193e45391caab9187a5def6f86f548.nc')

########### EXTRAINDO DADO DE SHAPEFILE #############
RS = ('C:/Users/feder/OneDrive/Área de Trabalho/geo/IC Chico/2025/RS_Municipios_2023.shp')
municipios = gpd.read_file(RS)
vale_taquari_mun = ['Arroio do Meio',
                    'Arvorezinha', 
                    'Boqueirão do Leão',
                    'Canudos do Vale',
                    'Capitão',
                    'Colinas',
                    'Coqueiro Baixo',
                    'Cruzeiro do Sul',
                    'Dois Lajeados',
                    'Encantado',
                    'Estrela',
                    'Forquetinha', 
                    'Ilópolis',
                    'Imigrante',
                    'Itapuca',
                    'Lajeado',
                    'Marques de Souza',
                    'Muçum',
                    'Nova Bréscia',
                    'Paverama',
                    'Poço das Antas',
                    'Pouso Novo',
                    'Progresso',
                    'Putinga',
                    'Relvado',
                    'Roca Sales',
                    'Santa Clara do Sul',
                    'São Valentim do Sul',
                    'Sério',
                    'Taquari',
                    'Travesseiro',
                    'Vespasiano Corrêa',
                    'Westfalia'
    ]
municipios_vale = municipios[municipios['NM_MUN'].isin(vale_taquari_mun)]
########### EXTRAINDO DADO DE SHAPEFILE #############

########### EXTRAINDO NETCDF ###########
ds_tp = xr.open_dataset(xr.backends.NetCDF4DataStore(tpnetcdf4))
ds_180 = ds_tp.assign_coords(longitude=(((ds_tp.longitude + 180) % 360) - 180)).sortby('longitude')
ds_180['tp']
ds_180['tp'].plot()
da_tp = ds_180['tp'] * 100000
da_tp[0,:,:].plot()
########### EXTRAINDO NETCDF ###########

tp = da_tp.sel(valid_time=da_tp.valid_time[0])
levels = np.linspace(200, 500, 13)

# Define a projeção do mapa
fig, ax = plt.subplots(figsize=(11, 6), subplot_kw={'projection': ccrs.PlateCarree()})

img = tp.plot(ax=ax, cmap='YlGnBu', transform=ccrs.PlateCarree(), levels=levels,
             extend='neither', cbar_kwargs={'label': 'Precipitação (mm)', 'shrink': 0.7}, zorder=1)

municipios.boundary.plot(ax=ax, edgecolor='black', linewidth=0.7,
                               transform=ccrs.PlateCarree(), zorder=2)
#ax.gridlines(draw_labels=dict(left=True, bottom=True, top=False, right=False))

for idx, row in municipios_vale.iterrows():
    centroid = row.geometry.centroid
    ax.text(centroid.x, centroid.y, row['NM_MUN'],
            fontsize=6, ha='center', va='center',
            transform=ccrs.PlateCarree(), zorder=3)


# Título
ax.set_title("Precipitação setembro 2023 - Vale do Taquari", fontsize=13)
ax.set_extent([-52.6, -51.5, -30, -28.5], crs=ccrs.PlateCarree())
# Eixos
#ax.gridlines
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

plt.tight_layout()
plt.show()