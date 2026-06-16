"""
Compare LandTrendr forest disturbance trends with AFV trajectories
per breeding biome.
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.signal import savgol_filter
import os

os.chdir(r"D:\1-PhD Research data\Pollination service\Estimating Avian Pollination Service in North America\validation\gee_export")
out_dir = r"D:\1-PhD Research data\Pollination service\Estimating Avian Pollination Service in North America\output"

# ---- 1. Load LandTrendr data ----
lt = pd.read_csv("LandTrendr_biome_forest_loss_1985_2017.csv")
lt['loss_km2'] = lt['sum'] / 1e6
lt['loss_ha'] = lt['sum'] / 1e4

# GEE biome names have underscore -> map to AFV convention
name_map = {'Western_Forest': 'Western Forest', 'Eastern_Forest': 'Eastern Forest', 'Aridlands': 'Aridlands'}
lt['biome_clean'] = lt['biome'].map(name_map)

print("=== Forest loss summary (km2) ===")
print(lt.groupby('biome_clean')['loss_km2'].describe())

# ---- 2. Load AFV data ----
cont = pd.read_csv(os.path.join(out_dir, "Continental level N and AFV trajectories.csv"))
biome_afv = pd.read_csv(os.path.join(out_dir, "Breeding.Biome level N and AFV trajectories.csv"))
sp = pd.read_csv(os.path.join(out_dir, "species level N and AFV trajectories.csv"))

BIOMES = ['Western Forest', 'Eastern Forest', 'Aridlands']

print(f"\nContinental years: {cont['year'].min():.0f}-{cont['year'].max():.0f}")
print(f"Biomes: {list(biome_afv['level'].unique())}")
print(f"Species: {sp['level'].nunique()}")

# ---- 3. Per-biome forest loss trends ----
print("\n=== Per-biome forest loss trends (1985-2017) ===")
biome_trends = {}
for biome in BIOMES:
    b = lt[lt['biome_clean'] == biome]
    slope, intercept, r, p, se = stats.linregress(b['year'], b['loss_km2'])
    first = b[b['year'] == b['year'].min()]['loss_km2'].values[0]
    last  = b[b['year'] == b['year'].max()]['loss_km2'].values[0]
    pct = (last - first) / first * 100
    biome_trends[biome] = {'slope': slope, 'p': p, 'r2': r**2, 'pct_change': pct}
    direction = 'UP' if slope > 0 else 'DOWN'
    sig = '**' if p < 0.01 else ('*' if p < 0.05 else '')
    print(f"  {biome:20s}: {direction:4s} {slope:+.1f} km2/yr, R2={r**2:.3f}, p={p:.4f} {sig}")

# ---- 4. AFV trends per biome ----
print("\n=== AFV trends per biome (1970-2017) ===")
afv_trends = {}
for biome in BIOMES:
    b = biome_afv[biome_afv['level'] == biome]
    slope, intercept, r, p, se = stats.linregress(b['year'], b['AFV_med'])
    first_val = b[b['year'] == 1970]['AFV_med'].values[0] / 1e9
    last_val = b[b['year'] == 2017]['AFV_med'].values[0] / 1e9
    pct = (last_val - first_val) / first_val * 100
    afv_trends[biome] = {'slope': slope, 'p': p, 'r2': r**2, 'pct_change': pct}
    direction = 'UP' if slope > 0 else 'DOWN'
    sig = '**' if p < 0.01 else ('*' if p < 0.05 else '')
    print(f"  {biome:20s}: {direction:4s} {slope/1e9:+.3f}B/yr, R2={r**2:.3f}, p={p:.4f} {sig}")

# ---- 5. Hypothesis evaluation ----
print("\n=== HYPOTHESIS EVALUATION ===")
print("Discussion claims: Western Forest early seral habitat decline -> Selasphorus AFV drop\n")
for biome in BIOMES:
    lt_t = biome_trends[biome]
    afv_t = afv_trends[biome]
    lt_dir = 'UP' if lt_t['slope'] > 0 else 'DOWN'
    afv_dir = 'UP' if afv_t['slope'] > 0 else 'DOWN'

    print(f"  {biome}:")
    print(f"    Forest loss:  {lt_dir} ({lt_t['slope']:+.0f} km2/yr, p={lt_t['p']:.4f})")
    print(f"    AFV trend:    {afv_dir} ({afv_t['slope']/1e9:+.3f}B/yr, p={afv_t['p']:.4f})")

    if biome == 'Western Forest':
        if lt_t['slope'] < 0:
            print(f"    >>> Forest disturbance DECLINING -> early seral habitat IS shrinking")
            print(f"    >>> CONSISTENT with 'fire suppression + logging decline -> less nectar' hypothesis")
        else:
            print(f"    >>> Forest disturbance INCREASING -> MORE early seral habitat")
            print(f"    >>> REFUTES 'early seral habitat loss' hypothesis")
    elif biome == 'Eastern Forest':
        print(f"    (Ruby-throated Hummingbird biome — most adaptable pollinator)")
    print()

# ---- 6. Spearman correlation ----
print("=== Spearman correlation: forest loss vs AFV (aligned 1985-2017) ===")
for biome in BIOMES:
    lt_b = lt[lt['biome_clean'] == biome][['year', 'loss_km2']]
    afv_b = biome_afv[biome_afv['level'] == biome][['year', 'AFV_med']]
    m = pd.merge(lt_b, afv_b, on='year', how='inner').dropna()
    r, p = stats.spearmanr(m['AFV_med'], m['loss_km2'])
    direction = 'positive' if r > 0 else 'negative'
    sig = ' *' if p < 0.05 else ''
    print(f"  {biome:20s}: rho={r:+.3f} ({direction}), p={p:.4f}, n={len(m)}{sig}")

# ---- 7. FIGURES ----
fig = plt.figure(figsize=(18, 12))
colors = {'Western Forest': '#D55E00', 'Eastern Forest': '#0072B2', 'Aridlands': '#009E73'}

# Panel A: Forest loss
ax = fig.add_subplot(2, 3, 1)
for biome in BIOMES:
    b = lt[lt['biome_clean'] == biome]
    ax.plot(b['year'], b['loss_km2'], color=colors[biome], lw=2, label=biome)
ax.set_xlabel('Year'); ax.set_ylabel('Forest loss (km2)')
ax.set_title('A: Forest fast-loss area (USFS LCMS 1985-2017)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panel B: AFV by biome
ax = fig.add_subplot(2, 3, 2)
for biome in BIOMES:
    b = biome_afv[biome_afv['level'] == biome]
    ax.fill_between(b['year'], b['AFV_lci']/1e9, b['AFV_uci']/1e9, alpha=0.12, color=colors[biome])
    ax.plot(b['year'], b['AFV_med']/1e9, color=colors[biome], lw=2, label=biome)
ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('Year'); ax.set_ylabel('AFV (Billions)')
ax.set_title('B: AFV trajectory by breeding biome (1970-2017)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panel C: % of 1985 baseline
ax = fig.add_subplot(2, 3, 3)
for biome in BIOMES:
    b = lt[lt['biome_clean'] == biome].sort_values('year')
    baseline = b[b['year'] == 1985]['loss_km2'].values[0]
    b = b.copy(); b['pct'] = b['loss_km2'] / baseline * 100
    ax.plot(b['year'], b['pct'], color=colors[biome], lw=2, label=biome)
ax.axhline(y=100, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('Year'); ax.set_ylabel('% of 1985 baseline')
ax.set_title('C: Forest loss relative to 1985')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Panels D-F: dual-axis
panels = [('D: Western Forest', 4), ('E: Eastern Forest', 5), ('F: Aridlands', 6)]
for title, pos in panels:
    biome = title.split(': ')[1]
    ax = fig.add_subplot(2, 3, pos)

    afv_b = biome_afv[biome_afv['level'] == biome]
    ax.fill_between(afv_b['year'], afv_b['AFV_lci']/1e9, afv_b['AFV_uci']/1e9,
                    alpha=0.1, color='#AA4466')
    ax.plot(afv_b['year'], afv_b['AFV_med']/1e9, color='#AA4466', lw=2)
    ax.set_ylabel('AFV (Billions)', color='#AA4466', fontsize=9)
    ax.tick_params(axis='y', labelcolor='#AA4466')

    lt_b = lt[lt['biome_clean'] == biome]
    ax2 = ax.twinx()
    ax2.bar(lt_b['year'], lt_b['loss_km2'], color='#0072B2', alpha=0.25, width=0.8)
    ax2.set_ylabel('Forest loss (km2)', color='#0072B2', fontsize=9)
    ax2.tick_params(axis='y', labelcolor='#0072B2')

    # Smoothed forest loss trend
    if len(lt_b) >= 7:
        smoothed = savgol_filter(lt_b['loss_km2'].values, window_length=7, polyorder=2)
        ax2.plot(lt_b['year'], smoothed, color='#0072B2', lw=2)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('Year')
    ax.set_xlim(1984, 2018)

    afv_s = afv_trends[biome]['slope'] / 1e9
    lt_s = biome_trends[biome]['slope']
    ax.text(0.98, 0.02,
            f'AFV: {"UP" if afv_s>0 else "DOWN"} {abs(afv_s):.2f}B/yr  '
            f'Loss: {"UP" if lt_s>0 else "DOWN"} {abs(lt_s):.0f}km2/yr',
            transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

plt.suptitle('LandTrendr Forest Disturbance vs Avian Pollinator Functional Value by Breeding Biome',
             fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('landtrendr_afv_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("\nSaved: landtrendr_afv_comparison.png")

# ---- 8. Selasphorus-specific ----
print("\n=== Selasphorus (Rufous + Allen's + Calliope) vs Western Forest ===")
sel_spp = ["Rufous Hummingbird", "Calliope Hummingbird", "Broad-tailed Hummingbird"]
sel = sp[sp['level'].isin(sel_spp)]
sel_afv = sel.groupby('year').agg(
    AFV_med=('AFV_med', 'sum'), AFV_lci=('AFV_lci', 'sum'), AFV_uci=('AFV_uci', 'sum')
).reset_index()

wf = lt[lt['biome_clean'] == 'Western Forest'][['year', 'loss_km2']]
m = pd.merge(sel_afv, wf, on='year', how='inner')

sel_slope, _, sel_r2, sel_p, _ = stats.linregress(m['year'], m['AFV_med'] / 1e9)
lt_slope, _, lt_r2, lt_p, _ = stats.linregress(m['year'], m['loss_km2'])
r_spear, p_spear = stats.spearmanr(m['AFV_med'], m['loss_km2'])

print(f"  Selasphorus AFV: {sel_slope:+.3f}B/yr, R2={sel_r2**2:.3f}, p={sel_p:.4f}")
print(f"  Western Forest loss: {lt_slope:+.1f}km2/yr, R2={lt_r2**2:.3f}, p={lt_p:.4f}")
print(f"  Spearman rho = {r_spear:+.3f}, p = {p_spear:.4f}, n = {len(m)}")

if sel_slope < 0 and lt_slope < 0:
    print("  >>> Both declining: CONSISTENT with early seral habitat loss hypothesis")
elif sel_slope < 0 and lt_slope > 0:
    print("  >>> AFV down, disturbance UP: REFUTES early seral habitat loss hypothesis")
else:
    print("  >>> Mixed pattern — check")

# Selasphorus plot
fig2, ax = plt.subplots(figsize=(12, 6))
ax.fill_between(m['year'], m['AFV_lci']/1e9, m['AFV_uci']/1e9, alpha=0.15, color='#D55E00')
ax.plot(m['year'], m['AFV_med']/1e9, color='#D55E00', lw=2.5,
        label="Selasphorus AFV (Rufous + Calliope + Broad-tailed)")
ax.set_ylabel('AFV (Billions)', color='#D55E00', fontsize=12)
ax.tick_params(axis='y', labelcolor='#D55E00')
ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)

ax2 = ax.twinx()
ax2.bar(m['year'], m['loss_km2'], color='#0072B2', alpha=0.25, width=0.8, label='Forest fast-loss')
if len(m) >= 7:
    smoothed = savgol_filter(m['loss_km2'].values, window_length=7, polyorder=2)
    ax2.plot(m['year'], smoothed, color='#0072B2', lw=2, label='Loss trend (Savitzky-Golay)')
ax2.set_ylabel('Forest loss (km2)', color='#0072B2', fontsize=12)
ax2.tick_params(axis='y', labelcolor='#0072B2')

lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)
ax.set_title(f'Selasphorus AFV vs Western Forest Disturbance (1985-2017)\n'
             f'Spearman rho = {r_spear:+.3f} (p = {p_spear:.3f})',
             fontsize=13)
ax.set_xlabel('Year'); ax.set_xlim(1984, 2018)
plt.tight_layout()
plt.savefig('selasphorus_vs_western_forest.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: selasphorus_vs_western_forest.png")

print("\n=== DONE ===")
