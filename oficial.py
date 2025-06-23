# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 12:19:13 2025

@author: Federico
"""
import cdsapi
import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import stats
import numpy as np


dataset = "reanalysis-era5-land-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": ["2m_temperature"],
    "year": [
        "1988", "1989", "1990",
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
        "2024"
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
    "area": [0, -80, -60, -45]
}

url= 'https://cds.climate.copernicus.eu/api'
key= 'a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a'
client = cdsapi.Client(url,key)
client.retrieve(dataset, request).download()

ds = xr.open_dataset('758478d747afa9f806b08dcf651dac8a.nc')
ds
t2m = ds['t2m'] - 273.15 

# Climatologia (1988–2023)
climatologia = t2m.sel(valid_time=slice("1988-01-01", "2023-12-01")).mean("valid_time")

# Plot
shp_path = "C:/Users/feder/OneDrive/Documentos/file_uh/UH.shp"

gdf = gpd.read_file(shp_path)

levels = np.linspace(-1.5, 1.5, 13)

# Criar a figura com 3 subplots lado a lado
for ano in range(1985, 2024):
    # Média anual do ano específico
    media_anual = t2m.sel(valid_time=slice(f"{ano}-01-01", f"{ano}-12-01")).mean("valid_time")

    # Anomalia anual
    anomalia_anual = media_anual - climatologia

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"projection": ccrs.PlateCarree()})

    contour = ax.contourf(
        anomalia_anual.longitude,
        anomalia_anual.latitude,
        anomalia_anual,
        cmap='coolwarm',
        extend='both',
        levels=levels,
        transform=ccrs.PlateCarree(),
    )

    gdf.boundary.plot(ax=ax, edgecolor="black", linewidth=1)
    #gdf_com.boundary.plot(ax=ax, edgecolor="blue", linewidth=8)
    ax.plot(-71, -13.9, marker="o", color="green", markersize=10, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.set_extent([-74, -70.5, -14.8, -10], crs=ccrs.PlateCarree())
    # Grade
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 8}
    gl.ylabel_style = {"size": 8}

    # Barra de cores
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05, aspect=30, shrink=0.7)
    cbar.set_label('(°C)', size=10)

    # Título com o ano
    ax.set_title(f"Anomalia de temperatura do ar - {ano}", fontsize=11)

    plt.tight_layout()
    plt.show()