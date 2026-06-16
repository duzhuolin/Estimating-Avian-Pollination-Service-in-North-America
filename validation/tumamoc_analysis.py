import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os

os.chdir(r"D:\1-PhD Research data\Pollination service\Estimating Avian Pollination Service in North America\validation")

# Load data (only density; species names are hardcoded below)
density = pd.read_csv("tumamoc_data/CsvFiles/SMDensity.csv")

print(f"SMDensity: {density.shape[0]} rows, {density['Code'].nunique()} species, {density['Year'].nunique()} years")
print(f"Years: {sorted(density['Year'].unique())}")

# Bird-pollinated species (based on known pollination syndromes)
bird_codes = {
    'FOSP': ('Ocotillo', 'Hummingbird'),
    'CAGI': ('Saguaro', 'Bat/Bird'),
    'FEWI': ('Barrel cactus', 'Hummingbird'),
    'ECFA': ('Hedgehog cactus', 'Hummingbird'),
    'CAER': ('Fairyduster', 'Hummingbird'),
    'HYEM': ('Desert lavender', 'Hummingbird'),
    'SACO': ('Chia', 'Hummingbird'),
    'PEPA': ("Parry's penstemon", 'Hummingbird'),
    'SILO': ('Longflower tubetongue', 'Hummingbird'),
    'DESC': ('Larkspur', 'Hummingbird'),
    'LYEX': ('Wolfberry', 'Hummingbird/Bee'),
    'LYBE': ('Wolfberry', 'Hummingbird/Bee'),
    'LYFR': ('Wolfberry', 'Hummingbird/Bee'),
    'HICO': ('Desert rosemallow', 'Hummingbird'),
}

# Filter density data
bird_density = density[density['Code'].isin(bird_codes.keys())].copy()
bird_density['CommonName'] = bird_density['Code'].map(lambda x: bird_codes[x][0])
bird_density['Pollinator'] = bird_density['Code'].map(lambda x: bird_codes[x][1])

print(f"\nBird-pollinated spp in density data: {bird_density['Code'].nunique()}")
print(f"Rows: {bird_density.shape[0]}")

# Check which species have data
print("\n--- Species present in density records ---")
for code in sorted(bird_codes.keys()):
    n = (bird_density['Code'] == code).sum()
    yrs = sorted(bird_density[bird_density['Code'] == code]['Year'].unique())
    print(f"  {code} ({bird_codes[code][0]}): {n} records, years: {yrs}")

# ---- Analysis 1: Annual total density (1970-2012) ----
bird_study = bird_density[bird_density['Year'] >= 1970].copy()

# Total count across all plots per species per year
annual = bird_study.groupby(['Code', 'CommonName', 'Year'])['Count'].agg(
    ['sum', 'mean', 'count']
).reset_index()

print(f"\n--- Study period (1970-2012) annual data ---")
print(f"Species with >= 3 years of data in study period:")
for code in sorted(bird_codes.keys()):
    sp_data = annual[annual['Code'] == code]
    if len(sp_data) >= 3:
        print(f"  {code} ({bird_codes[code][0]}): {len(sp_data)} years")

# ---- Analysis 2: Trend per species ----
print("\n--- Per-species density trends (1970-2012) ---")
trends = []
for code in sorted(bird_codes.keys()):
    sp_data = annual[annual['Code'] == code].dropna()
    if len(sp_data) >= 3:
        slope, intercept, r, p, se = stats.linregress(sp_data['Year'], sp_data['sum'])
        first_val = sp_data[sp_data['Year'] == sp_data['Year'].min()]['sum'].values[0]
        last_val  = sp_data[sp_data['Year'] == sp_data['Year'].max()]['sum'].values[0]
        trends.append({
            'Code': code, 'CommonName': bird_codes[code][0],
            'slope': slope, 'p_value': p, 'r_squared': r**2,
            'first_total': first_val, 'last_total': last_val,
            'n_years': len(sp_data),
            'direction': 'Increase' if slope > 0 else 'Decrease'
        })
        sig = '**' if p < 0.05 else ('*' if p < 0.1 else '')
        print(f"  {code:5s} {bird_codes[code][0]:25s} slope={slope:+.3f}/yr p={p:.3f} R²={r**2:.3f} {sig}")

# ---- Analysis 3: Aggregate total ----
total_annual = bird_study.groupby('Year')['Count'].sum().reset_index()
total_annual['n_species'] = bird_study.groupby('Year')['Code'].nunique().values

if len(total_annual) >= 3:
    slope, intercept, r, p, se = stats.linregress(total_annual['Year'], total_annual['Count'])
    print(f"\n--- Aggregate trend (all bird-pollinated species, 1970-2012) ---")
    print(f"  Total individuals counted: slope={slope:+.3f}/yr, p={p:.3f}, R²={r**2:.3f}")
    first_tot = total_annual[total_annual['Year'] == total_annual['Year'].min()]['Count'].values[0]
    last_tot  = total_annual[total_annual['Year'] == total_annual['Year'].max()]['Count'].values[0]
    print(f"  {int(total_annual['Year'].min())}: {first_tot:.0f} -> {int(total_annual['Year'].max())}: {last_tot:.0f}")

# ---- Plot ----
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# Panel A: Per-species total counts
ax = axes[0]
for code in sorted(bird_codes.keys()):
    sp_data = annual[annual['Code'] == code].dropna()
    if len(sp_data) >= 3:
        ax.plot(sp_data['Year'], sp_data['sum'], 'o-', label=f"{code} ({bird_codes[code][0]})",
                lw=1.5, markersize=4)
ax.set_xlabel('Year')
ax.set_ylabel('Total individuals (all plots)')
ax.set_title('Hummingbird-pollinated perennials — Tumamoc Hill (1970–2012)')
ax.legend(fontsize=7, ncol=2, loc='upper left')
ax.grid(True, alpha=0.3)

# Panel B: Aggregate
ax = axes[1]
ax.plot(total_annual['Year'], total_annual['Count'], 'o-', color='#31688EFF', lw=2, markersize=8)
ax.set_xlabel('Year')
ax.set_ylabel('Total individuals (all species)')
ax.set_title('Aggregate density — all hummingbird-pollinated species')
ax.grid(True, alpha=0.3)
if len(total_annual) >= 3:
    x_range = np.linspace(total_annual['Year'].min(), total_annual['Year'].max(), 100)
    ax.plot(x_range, intercept + slope * x_range, 'r--', lw=1, alpha=0.5,
            label=f'slope={slope:+.2f}/yr, p={p:.3f}')
    ax.legend()

plt.tight_layout()
plt.savefig('tumamoc_bird_pollinated_density_trends.png', dpi=150)
plt.close()
print("\nPlot saved to validation/tumamoc_bird_pollinated_density_trends.png")
print("\n=== DONE ===")
