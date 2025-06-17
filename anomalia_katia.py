import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd


# ==============================================
# 1. Configuração inicial
# ==============================================
# Ajuste o caminho para seu arquivo ERA5
path_era5 = "C:/Users/Federico/Downloads/data_stream-moda_stepType-avgua.nc"  # Substitua pelo seu caminho!

# Carregar dados
ds = xr.open_dataset(path_era5)
t2m_kelvin = ds['t2m']  # Temperatura em Kelvin
tp = ds['tp']           # Precipitação em metros

# Converter Kelvin para Celsius
t2m_celsius = t2m_kelvin - 273.15

# Definir períodos decadais
decadas = {
    '1988-2000': slice('1988-01-01', '2000-12-31'),
    '2001-2010': slice('2001-01-01', '2010-12-31'),
    '2011-2023': slice('2011-01-01', '2023-12-31')
}

# Referência climatológica (média 1988-2023 em °C)
ref_t2m = t2m_celsius.sel(valid_time=slice('1988-01-01', '2023-12-31')).mean('valid_time')

# ==============================================
# 2. Calcular anomalias para t2m (em °C)
# ==============================================
anomalias_t2m = {}
for nome, periodo in decadas.items():
    media_decada = t2m_celsius.sel(valid_time=periodo).mean('valid_time')
    anomalias_t2m[nome] = media_decada - ref_t2m  # Anomalia em °C

# ==============================================
# 3. Plotar anomalias de t2m (3 subplots lado a lado)
# ==============================================

path_shp = 'C:/Users/Federico/Downloads/Federico/bacia - Copia.shp'
bacia = gpd.read_file(path_shp)

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 6))

# Função para configurar mapas
def configure_map(ax, anom_data, title):
    # Plotar shapefile da bacia
    bacia.plot(ax=ax, color='none', edgecolor='black', linewidth=1.5)
    
    # Plotar anomalias
    anom_data.plot(
        ax=ax,
        cmap='coolwarm',
        vmin=-2, vmax=2,
        add_colorbar=False
    )
    ax.set_title(title, fontsize=12)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

# Plotar cada década
for i, (nome, anom) in enumerate(anomalias_t2m.items()):
    ax = axes[i]
    configure_map(ax, anom, f'T2M Anomalia ({nome}) [°C]')

# Adicionar barra de cores única
cbar = fig.colorbar(
    plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=-2, vmax=2)),
    ax=axes,
    orientation='horizontal',
    pad=0.1,
    shrink=0.8,
    aspect=40,
    label='Anomalia de Temperatura (°C)'
)

plt.suptitle('Anomalias Decadais de Temperatura (ERA5, 1988-2023)', fontsize=14, y=1.05)
plt.show()
