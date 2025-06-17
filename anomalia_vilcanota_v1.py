# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 19:47:35 2025

@author: federico
"""


import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import stats

# === 1. Abrir o NetCDF e converter K para °C ===
ds = xr.open_dataset("C:/Users/feder/OneDrive/Documentos/data_stream-moda_stepType-avgua.nc")
t2m = ds['t2m'] - 273.15  # Kelvin para Celsius

# === 2. Calcular a climatologia (1982-2023) ===
clim = t2m.sel(valid_time=slice("1982", "2023")).groupby("valid_time.month").mean("valid_time")
clim

print(t2m.dims)

# === 3. Remover tendência linear ===
# Regressão linear em cada ponto da grade
def remove_trend(da):
    time_num = np.arange(da.sizes['valid_time'])

    def regress(y):
        slope, intercept, *_ = stats.linregress(time_num, y)
        return y - (slope * time_num + intercept)

    return xr.apply_ufunc(
        regress, da,
        input_core_dims=[['valid_time']],
        output_core_dims=[['valid_time']],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[da.dtype]
    )

t2m_detrended = remove_trend(t2m)


# === 4. Calcular as anomalias decadais ===
decadas = {
    "1988-2000": ("1988", "2000"),
    "2001-2010": ("2001", "2010"),
    "2011-2023": ("2011", "2023"),
}

anomalias = {}
for nome, (inicio, fim) in decadas.items():
    dados = t2m_detrended.sel(valid_time=slice(inicio, fim))
    media = dados.groupby("valid_time.month").mean("valid_time")
    anomalia = media - clim
    anomalias[nome] = anomalia

# === 5. Carregar o shapefile ===
shapefile_path = r"C:/Users/feder/OneDrive/Documentos/Federico/bacia - Copia.shp"
gdf = gpd.read_file(shapefile_path)

# === 6. Plotar com uma barra de cores única ===
fig, axs = plt.subplots(1, 3, figsize=(18, 6), subplot_kw={'projection': ccrs.PlateCarree()})
vmin, vmax = -2.5, 2.5  # defina um intervalo comum para a barra de cores

for ax, (nome, anom) in zip(axs, anomalias.items()):
    anom_mean = anom.mean('month')  # média anual da anomalia decadal
    im = anom_mean.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap="coolwarm",
        vmin=vmin,
        vmax=vmax,
        cbar_kwargs={'label': 'Anomalia (°C)', 'shrink': 0.7}
    )
    gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=1)
    ax.coastlines()
    ax.set_title(f"Anomalia: {nome}")
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray', alpha=0.3)

plt.tight_layout()
plt.show()
