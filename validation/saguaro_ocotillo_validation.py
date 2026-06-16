"""
Compare Tumamoc Hill saguaro/ocotillo density trends with
species-specific AFV from the pollination model.

Saguaro pollinators: White-winged Dove + Cactus Wren
Ocotillo pollinators: Aridlands hummingbirds (Anna's, Black-chinned, Costa's)
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# ---- 1. Load species-level AFV data ----
afv_path = r"D:\1-PhD Research data\Pollination service\Estimating Avian Pollination Service in North America\output\species level N and AFV trajectories.csv"
afv = pd.read_csv(afv_path)

# Check available species
print("Available species:", sorted(afv['level'].unique()))

# ---- 2. Saguaro analysis: White-winged Dove + Cactus Wren ----
wwdo = afv[afv['level'] == 'White-winged Dove'].copy()
cawr = afv[afv['level'] == 'Cactus Wren'].copy()

# Sum AFV across both pollinators
saguaro_pollinators = pd.merge(
    wwdo[['year', 'AFV_med', 'AFV_lci', 'AFV_uci']],
    cawr[['year', 'AFV_med', 'AFV_lci', 'AFV_uci']],
    on='year', suffixes=('_wwdo', '_cawr')
)
saguaro_pollinators['AFV_med'] = saguaro_pollinators['AFV_med_wwdo'] + saguaro_pollinators['AFV_med_cawr']
saguaro_pollinators['AFV_lci'] = saguaro_pollinators['AFV_lci_wwdo'] + saguaro_pollinators['AFV_lci_cawr']
saguaro_pollinators['AFV_uci'] = saguaro_pollinators['AFV_uci_wwdo'] + saguaro_pollinators['AFV_uci_cawr']

# Saguaro density data from Tumamoc Hill (1975-2012)
saguaro_density = pd.DataFrame({
    'Year': [1975, 1978, 1984, 1985, 1993, 2001, 2010, 2012],
    'mean_per_plot': [3.0, 2.3, 1.7, 1.0, 2.0, 2.6, 2.4, 2.0],
    'total_count': [6, 7, 5, 1, 2, 13, 12, 2]
})

# ---- 3. Ocotillo analysis: Aridlands hummingbirds ----
# Identify Aridlands-breeding hummingbirds from the species list
# Anna's, Black-chinned, Costa's are the core Aridlands hummingbirds
arid_hummers = ['Anna\'s Hummingbird', "Allen's Hummingbird", "Costa's Hummingbird"]

ocotillo_pollinators = afv[afv['level'].isin(arid_hummers)]
ocotillo_afv = ocotillo_pollinators.groupby('year').agg(
    AFV_med = ('AFV_med', 'sum'),
    AFV_lci = ('AFV_lci', 'sum'),
    AFV_uci = ('AFV_uci', 'sum')
).reset_index()

# Ocotillo density data from Tumamoc Hill (1974-2012)
ocotillo_density = pd.DataFrame({
    'Year': [1974, 1975, 1978, 1984, 1985, 2001, 2010, 2012],
    'mean_per_plot': [1.0, 1.0, 1.67, 1.0, 2.0, 1.33, 1.0, 2.5],
    'total_count': [1, 1, 5, 1, 2, 4, 3, 5]
})

# ---- 4. Print comparison ----
print("\n=== SAGUARO: White-winged Dove + Cactus Wren AFV vs Saguaro density ===")
for _, row in saguaro_density.iterrows():
    yr = row['Year']
    afv_row = saguaro_pollinators[saguaro_pollinators['year'] == yr]
    if len(afv_row) > 0:
        afv_val = afv_row['AFV_med'].values[0] / 1e9
        print(f"  {yr}: AFV={afv_val:.1f}B, saguaro density={row['mean_per_plot']:.1f}/plot (n={int(row['total_count'])})")

print("\n=== OCOTILLO: Aridlands hummingbird AFV vs Ocotillo density ===")
for _, row in ocotillo_density.iterrows():
    yr = row['Year']
    afv_row = ocotillo_afv[ocotillo_afv['year'] == yr]
    if len(afv_row) > 0:
        afv_val = afv_row['AFV_med'].values[0] / 1e9
        print(f"  {yr}: AFV={afv_val:.1f}B, ocotillo density={row['mean_per_plot']:.1f}/plot (n={int(row['total_count'])})")

# ---- 5. Plot: dual-panel ----
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Panel A: Saguaro
ax = axes[0]
ax.fill_between(saguaro_pollinators['year'],
                saguaro_pollinators['AFV_lci'] / 1e9,
                saguaro_pollinators['AFV_uci'] / 1e9,
                alpha=0.2, color='#D55E00')
ax.plot(saguaro_pollinators['year'], saguaro_pollinators['AFV_med'] / 1e9,
        color='#D55E00', lw=2, label='WW Dove + Cactus Wren AFV')
ax.set_ylabel('AFV (Billions)', color='#D55E00', fontsize=12)
ax.tick_params(axis='y', labelcolor='#D55E00')

ax2 = ax.twinx()
ax2.plot(saguaro_density['Year'], saguaro_density['mean_per_plot'],
         'o-', color='#009E73', lw=2, markersize=10, markerfacecolor='white',
         markeredgewidth=2, label='Saguaro density')
ax2.set_ylabel('Saguaro mean density (indiv/plot)', color='#009E73', fontsize=12)
ax2.tick_params(axis='y', labelcolor='#009E73')

ax.set_title('Saguaro pollinators (WW Dove + Cactus Wren) AFV vs. Saguaro density at Tumamoc Hill',
             fontsize=13)
ax.set_xlim(1970, 2017)

# Add lines from both
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)

# Panel B: Ocotillo
ax = axes[1]
ax.fill_between(ocotillo_afv['year'],
                ocotillo_afv['AFV_lci'] / 1e9,
                ocotillo_afv['AFV_uci'] / 1e9,
                alpha=0.2, color='#0072B2')
ax.plot(ocotillo_afv['year'], ocotillo_afv['AFV_med'] / 1e9,
        color='#0072B2', lw=2, label='Aridlands hummingbird AFV')
ax.set_ylabel('AFV (Billions)', color='#0072B2', fontsize=12)
ax.tick_params(axis='y', labelcolor='#0072B2')

ax3 = ax.twinx()
ax3.plot(ocotillo_density['Year'], ocotillo_density['mean_per_plot'],
         'o-', color='#CC79A7', lw=2, markersize=10, markerfacecolor='white',
         markeredgewidth=2, label='Ocotillo density')
ax3.set_ylabel('Ocotillo mean density (indiv/plot)', color='#CC79A7', fontsize=12)
ax3.tick_params(axis='y', labelcolor='#CC79A7')

ax.set_title('Aridlands hummingbirds (Anna\'s + Allen\'s + Costa\'s) AFV vs. Ocotillo density at Tumamoc Hill',
             fontsize=13)
ax.set_xlim(1970, 2017)

lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax3.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)

plt.tight_layout()
plt.savefig('saguaro_ocotillo_afv_validation.png', dpi=150, bbox_inches='tight')
plt.close()

# ---- 6. Correlation ----
print("\n--- Spearman correlation ---")
# Merge saguaro
sg = pd.merge(saguaro_pollinators, saguaro_density, left_on='year', right_on='Year', how='inner')
if len(sg) >= 4:
    r, p = stats.spearmanr(sg['AFV_med'], sg['mean_per_plot'])
    print(f"Saguaro: rho={r:.2f}, p={p:.3f}, n={len(sg)} pairs")

oc = pd.merge(ocotillo_afv, ocotillo_density, left_on='year', right_on='Year', how='inner')
if len(oc) >= 4:
    r, p = stats.spearmanr(oc['AFV_med'], oc['mean_per_plot'])
    print(f"Ocotillo: rho={r:.2f}, p={p:.3f}, n={len(oc)} pairs")

print("\nPlot saved to validation/saguaro_ocotillo_afv_validation.png")
print("=== DONE ===")
