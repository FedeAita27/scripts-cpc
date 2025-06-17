import cdsapi
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import geopandas as gpd
from statsmodels.tsa.seasonal import STL
import zipfile
import os

# ========= ETAPA 1: Download ========= #
client = cdsapi.Client(url='https://cds.climate.copernicus.eu/api',
                       key='a9d91cd4-53fb-4494-a7b6-e4e7e67e0e2a')

request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": ["2m_temperature", "total_precipitation"],
    "year": [str(y) for y in range(1988, 2025)],
    "month": [f"{m:02d}" for m in range(1, 13)],
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [-13, -72, -15, -70]
}

client.retrieve("reanalysis-era5-single-levels-monthly-means", request).download()

zip_path = "baa8ab3c3da125b19c734e1df80c65ad.zip"
extract_path = "./dados_era5"
os.makedirs(extract_path, exist_ok=True)

with zipfile.ZipFile(zip_path, 'r') as zip_ref:
    zip_ref.extractall(extract_path)

# Lista os arquivos .nc extraídos
arquivos_extraidos = [f for f in os.listdir(extract_path) if f.endswith(".nc")]
print("Arquivos extraídos:", arquivos_extraidos)

# ========= ETAPA 2: Identificar os arquivos pelo conteúdo ========= #
arquivo_t2m = None
arquivo_tp = None

for file in arquivos_extraidos:
    file_path = os.path.join(extract_path, file)
    try:
        ds_temp = xr.open_dataset(file_path)
        variaveis = list(ds_temp.data_vars)
        print(f"{file} contém variáveis: {variaveis}")
        if "t2m" in variaveis:
            arquivo_t2m = file_path
        elif "tp" in variaveis:
            arquivo_tp = file_path
        ds_temp.close()
    except Exception as e:
        print(f"Erro ao abrir {file}: {e}")

# Verificação final
if arquivo_t2m is None:
    raise ValueError("Não foi possível identificar o arquivo de temperatura (t2m).")
if arquivo_tp is None:
    raise ValueError("Não foi possível identificar o arquivo de precipitação (tp).")

# ========= ETAPA 3: Abrir os arquivos corretamente ========= #
ds_t2m = xr.open_dataset(arquivo_t2m)
ds_tp = xr.open_dataset(arquivo_tp)

print("Arquivos abertos com sucesso!")
# ========= ETAPA 3: Preprocessamento ========= #
# Renomear variáveis
t2m = ds_t2m['t2m'] - 273.15  # Convertendo para °C
tp = ds_tp['tp'] * 1000       # Convertendo para mm

# Usar 'valid_time' como dimensão temporal
t2m['time'] = t2m['valid_time']
tp['time'] = tp['valid_time']

# ========= ETAPA 4: Remoção de tendência e sazonalidade ========= #
def deseasonalize_detrend(data):
    """Remove tendência e sazonalidade com STL"""
    
    time = data['valid_time'].values  # pegamos a dimensão temporal fora
    
    def apply_stl(ts_1d):
        # Cria índice de tempo
        ts = pd.Series(ts_1d, index=pd.date_range(start=str(time[0]), periods=len(time), freq="MS"))
        ts_interp = ts.interpolate()
        stl = STL(ts_interp, period=12, robust=True)
        result = stl.fit()
        return result.resid.values
    
    return xr.apply_ufunc(
        apply_stl,
        data,
        input_core_dims=[["valid_time"]],
        output_core_dims=[["valid_time"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float]
    )

# Aplicando corretamente:
t2m_deseas = deseasonalize_detrend(ds_t2m['t2m'])
tp_deseas = deseasonalize_detrend(ds_tp['tp'])

# ========= ETAPA 5: Cálculo da climatologia e anomalias ========= #
def compute_decadal_anomalies(data, full_period=('1988-01-01', '2023-12-31')):
    climatology = data.sel(valid_time=slice(*full_period)).groupby("valid_time.month").mean("valid_time")
    
    def anomaly_for_decade(start, end):
        decadal_data = data.sel(valid_time=slice(start, end))
        decadal_monthly = decadal_data.groupby("valid_time.month").mean("valid_time")
        return decadal_monthly - climatology

    anom_1988_2000 = anomaly_for_decade("1988-01-01", "2000-12-31")
    anom_2001_2010 = anomaly_for_decade("2001-01-01", "2010-12-31")
    anom_2011_2023 = anomaly_for_decade("2011-01-01", "2023-12-31")

    return anom_1988_2000, anom_2001_2010, anom_2011_2023

anom_t2m_88_00, anom_t2m_01_10, anom_t2m_11_23 = compute_decadal_anomalies(t2m_deseas)
anom_tp_88_00, anom_tp_01_10, anom_tp_11_23 = compute_decadal_anomalies(tp_deseas)

# ========= ETAPA 6: Plotar com shapefile ========= #
# Ler o shapefile
shp_path = "C:/Users/feder/OneDrive/Documentos/bacia - CopiaUTM.shp"
gdf = gpd.read_file(shp_path)

# Plotar exemplo para temperatura 2011–2023 (janeiro)
plt.figure(figsize=(10, 8))
ax = plt.subplot(1, 1, 1)
anom_t2m_11_23.sel(month=1).plot(ax=ax, cmap="RdBu_r", cbar_kwargs={"label": "Anomalia T2M (°C)"})
gdf.boundary.plot(ax=ax, edgecolor="black")
plt.title("Anomalia Decadal T2M - Janeiro (2011–2023)")
plt.show()
