"""
Intensity Analysis with Regional Decomposition
"""

import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
import seaborn as sns

# -------------------------------------------------------------------
# 0. CONFIGURATION
# -------------------------------------------------------------------

matrix_folder = r"africa_5r_my_matrix"  #please change this to the folder where all the matrices are downloaded using "gee_code/01_matrix_generation.js" code
output_folder = r"IA_regional_decomposition"
os.makedirs(output_folder, exist_ok=True)

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'

LC_CLASSES = {
    1: 'CRP', 2: 'FST', 3: 'SHR', 4: 'GRS', 5: 'TUD',
    6: 'WET', 7: 'IMP', 8: 'BAL', 9: 'WTR', 10: 'PSI'
}

LC_COLORS = {
    'CRP': '#E6B800', 'FST': '#228B22', 'SHR': '#808000', 'GRS': '#32CD32',
    'TUD': '#A9A9A9', 'WET': '#00CED1', 'IMP': '#DC143C', 'BAL': '#8B4513',
    'WTR': '#4169E1', 'PSI': '#1E90FF'
}

REGION_COLORS = {
    'EAF': '#e41a1c', 'MED': '#377eb8', 'SAF': '#4daf4a',
    'SAH': '#984ea3', 'WAF': '#ff7f00'
}
REGION_ORDER = ['EAF', 'MED', 'SAF', 'SAH', 'WAF']

def save_table(df, folder, name):
    if df.empty: return
    df.to_csv(os.path.join(folder, f"{name}.csv"), index=False)

# -------------------------------------------------------------------
# 1. LOAD DATA
# -------------------------------------------------------------------

csv_files = [f for f in os.listdir(matrix_folder) if f.lower().endswith(".csv")]
dfs = []
for fname in csv_files:
    df = pd.read_csv(os.path.join(matrix_folder, fname))
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_')
    needed = ["region", "year_initial", "year_final", "interval_duration_yrs",
              "from_class", "to_class", "area_km2"]
    if all(c in df.columns for c in needed):
        dfs.append(df[needed])

all_df = pd.concat(dfs, ignore_index=True)

# Create AFRICA aggregate
africa_agg = all_df.groupby(['year_initial', 'year_final', 'interval_duration_yrs',
                             'from_class', 'to_class'])['area_km2'].sum().reset_index()
africa_agg['region'] = 'AFRICA'
all_df = pd.concat([all_df, africa_agg], ignore_index=True)

# -------------------------------------------------------------------
# 2. MATRIX BUILDER
# -------------------------------------------------------------------

def build_matrices(df_site):
    intervals = df_site[["year_initial", "year_final", "interval_duration_yrs"]].drop_duplicates().sort_values("year_initial")
    matrices = {}
    for _, r in intervals.iterrows():
        yi, yf, dur = int(r["year_initial"]), int(r["year_final"]), float(r["interval_duration_yrs"])
        sub = df_site[(df_site["year_initial"] == yi) & (df_site["year_final"] == yf)]
        mat = sub.pivot_table(index="from_class", columns="to_class",
                              values="area_km2", aggfunc="sum", fill_value=0.0)
        gc = sorted(LC_CLASSES.keys())
        mat = mat.reindex(index=gc, columns=gc, fill_value=0.0)
        matrices[(yi, yf)] = {"matrix": mat, "duration": dur}
    return matrices

def get_region_matrix(all_df, region, yi, yf):
    """Extract transition matrix for a specific region and interval"""
    sub = all_df[(all_df['region'] == region) &
                 (all_df['year_initial'] == yi) &
                 (all_df['year_final'] == yf)]
    mat = sub.pivot_table(index="from_class", columns="to_class",
                          values="area_km2", aggfunc="sum", fill_value=0.0)
    gc = sorted(LC_CLASSES.keys())
    return mat.reindex(index=gc, columns=gc, fill_value=0.0)

# -------------------------------------------------------------------
# 3. INTERVAL LEVEL with Regional Decomposition
# -------------------------------------------------------------------

def compute_interval_metrics(matrices):
    rows = []
    for (yi, yf), info in matrices.items():
        mat, dur = info["matrix"], info["duration"]
        total = mat.values.sum()
        change = total - np.trace(mat.values)
        intensity = (change / total) / dur * 100 if total > 0 else 0
        change_pct = (change / total) * 100 if total > 0 else 0
        rows.append({
            "interval": f"{yi}-{yf}", "year_initial": yi, "year_final": yf,
            "duration": dur, "total_area": total, "change_area": change,
            "intensity": intensity, "change_pct": change_pct
        })
    df = pd.DataFrame(rows).sort_values("year_initial", ascending=False)
    U = np.nan
    if not df.empty and df["duration"].sum() > 0:
        U = (df["change_area"].sum() / df["total_area"].mean()) / df["duration"].sum() * 100
    return df, U

def compute_interval_with_regional(matrices, all_df_ref):
    """Compute interval metrics with regional breakdown for CSV export"""
    regions = [r for r in REGION_ORDER if r in all_df_ref['region'].unique()]

    rows = []
    for (yi, yf), info in matrices.items():
        mat, dur = info["matrix"], info["duration"]
        total = mat.values.sum()
        change = total - np.trace(mat.values)
        intensity = (change / total) / dur * 100 if total > 0 else 0
        change_pct = (change / total) * 100 if total > 0 else 0

        row = {
            "interval": f"{yi}-{yf}",
            "year_initial": yi,
            "year_final": yf,
            "duration_yrs": dur,
            "total_area_km2": round(total, 2),
            "change_area_km2": round(change, 2),
            "intensity_pct_per_yr": round(intensity, 4),
            "change_pct_of_map": round(change_pct, 4)
        }

        # Add regional breakdown
        for reg in regions:
            reg_mat = get_region_matrix(all_df_ref, reg, yi, yf)
            reg_change = reg_mat.values.sum() - np.trace(reg_mat.values)
            reg_contrib_pct = (reg_change / change * 100) if change > 0 else 0
            reg_total = reg_mat.values.sum()
            reg_intensity = (reg_change / reg_total) / dur * 100 if reg_total > 0 else 0

            row[f"{reg}_change_km2"] = round(reg_change, 2)
            row[f"{reg}_contrib_pct"] = round(reg_contrib_pct, 2)
            row[f"{reg}_intensity"] = round(reg_intensity, 4)

        rows.append(row)

    df = pd.DataFrame(rows).sort_values("year_initial")

    # Calculate uniform intensity
    U = np.nan
    if not df.empty and df["duration_yrs"].sum() > 0:
        U = (df["change_area_km2"].sum() / df["total_area_km2"].mean()) / df["duration_yrs"].sum() * 100

    return df, U

def plot_interval_level(df_int, U, site_name, all_df_ref=None):
    """
    3-panel interval plot:
    - Panel 1: Interval Change Area (% of map) - SIZE
    - Panel 2: Annual Intensity (% per year) - SPEED
    - Panel 3: Regional Contribution (% of change) - WHERE [AFRICA only]
    """
    if df_int.empty: return

    # Sort ascending for plotting (oldest at bottom)
    df_int = df_int.sort_values('year_initial', ascending=True)

    intervals = df_int["interval"].tolist()
    y_pos = np.arange(len(intervals))
    is_africa = (site_name == "AFRICA") and (all_df_ref is not None)

    if is_africa:
        fig, axes = plt.subplots(1, 3, figsize=(16, 9), sharey=True,
                                 gridspec_kw={'wspace': 0, 'width_ratios': [1, 1, 0.7]})
        ax_size, ax_speed, ax_where = axes
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 9), sharey=True,
                                 gridspec_kw={'wspace': 0})
        ax_size, ax_speed = axes

    # --- PANEL 1: SIZE (Change Area %) ---
    ax_size.barh(y_pos, df_int["change_pct"], color='#d9d9d9', edgecolor='black', height=0.7)
    ax_size.invert_xaxis()
    ax_size.set_xlabel("Interval Change Area\n(% of map)", fontweight='bold', fontsize=12)
    ax_size.set_yticks(y_pos)
    ax_size.set_yticklabels(intervals, fontweight='bold', fontsize=12)
    ax_size.grid(axis='x', linestyle='--', alpha=0.5)
    ax_size.tick_params(axis='y', length=0)
    for t in ax_size.get_xticklabels(): t.set_fontweight('bold')

    # --- PANEL 2: SPEED (Annual Intensity) ---
    ax_speed.barh(y_pos, df_int["intensity"], color='#696969', edgecolor='black', height=0.7)
    if not np.isnan(U):
        ax_speed.axvline(U, color='red', linestyle='--', linewidth=2, label=f'Uniform: {U:.2f}%')
        ax_speed.text(U, len(intervals) - 0.3, f'U = {U:.2f}%', ha='center', va='bottom',
                      fontsize=11, color='red', fontweight='bold')
    ax_speed.set_xlabel("Annual Change Intensity\n(% per year)", fontweight='bold', fontsize=12)
    ax_speed.grid(axis='x', linestyle='--', alpha=0.5)
    ax_speed.tick_params(axis='y', length=0)
    for t in ax_speed.get_xticklabels(): t.set_fontweight('bold')

    # Spine styling
    ax_size.spines['right'].set_linewidth(1.5)
    ax_speed.spines['left'].set_visible(False)

    # --- PANEL 3: WHERE (Regional Contribution) [AFRICA only] ---
    if is_africa:
        ax_speed.spines['right'].set_linewidth(1.5)
        ax_where.spines['left'].set_visible(False)
        ax_where.tick_params(axis='y', length=0)

        # Calculate regional contributions
        left_offset = np.zeros(len(intervals))
        regions = [r for r in REGION_ORDER if r in all_df_ref['region'].unique()]

        for reg in regions:
            contrib = []
            for interval in intervals:
                yi, yf = map(int, interval.split('-'))
                africa_change = df_int[df_int['interval'] == interval]['change_area'].values[0]

                reg_mat = get_region_matrix(all_df_ref, reg, yi, yf)
                reg_change = reg_mat.values.sum() - np.trace(reg_mat.values)
                pct = (reg_change / africa_change * 100) if africa_change > 0 else 0
                contrib.append(pct)

            ax_where.barh(y_pos, contrib, left=left_offset, height=0.7,
                         color=REGION_COLORS[reg], edgecolor='white', linewidth=0.5, label=reg)
            left_offset += np.array(contrib)

        ax_where.set_xlim(0, 100)
        ax_where.set_xlabel("Regional Contribution\n(% of total change)", fontweight='bold', fontsize=12)
        ax_where.grid(axis='x', linestyle=':', alpha=0.5)
        for t in ax_where.get_xticklabels(): t.set_fontweight('bold')

        # Legend below
        handles, labels = ax_where.get_legend_handles_labels()
        ax_where.legend(handles[::-1], labels[::-1], title='Region', loc='upper center',
                       bbox_to_anchor=(0.5, -0.10), ncol=len(regions), frameon=False, fontsize=10)

    plt.suptitle(f"Interval Level Intensity Analysis: {site_name}", fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90, bottom=0.18 if is_africa else 0.12)
    plt.savefig(os.path.join(output_folder, f"L1_Interval_{site_name}.png"), bbox_inches='tight')
    plt.close()

# -------------------------------------------------------------------
# 4. CATEGORY LEVEL with Regional Decomposition
#    MODIFIED: Panels 1 & 2 merged at zero (butterfly chart style)
# -------------------------------------------------------------------

def human_format(x, pos):
    if x == 0: return ''  # Don't show zero here, we'll add it manually
    if abs(x) >= 1e6: return f'{x/1e6:.1f}M'
    if abs(x) >= 1e3: return f'{x/1e3:.0f}k'
    return f'{x:.0f}'

def human_format_abs(x, pos):
    """Format showing absolute values (for left side of butterfly chart)"""
    if x == 0: return ''
    x = abs(x)
    if x >= 1e6: return f'{x/1e6:.1f}M'
    if x >= 1e3: return f'{x/1e3:.0f}k'
    return f'{x:.0f}'

def analyze_category_level(matrices, site_name, all_df_ref=None, grayscale_bars=False):
    """
    MODIFIED: Butterfly chart style - Area (LEFT) | Intensity (RIGHT) meeting at ZERO
    - Panel 1 (Area) and Panel 2 (Intensity) merged with shared zero axis
    - Panel 3: Regional Contribution to Gain/Loss - WHERE [AFRICA only]

    Parameters:
    - grayscale_bars: If True, panels 1 & 2 use gray colors instead of LC colors
    """
    print(f"   [Level 2] Category Analysis: {site_name}" + (" (grayscale)" if grayscale_bars else ""))

    lc_ids = sorted(LC_CLASSES.keys())
    classes = [LC_CLASSES[i] for i in lc_ids]
    is_africa = (site_name == "AFRICA") and (all_df_ref is not None)
    regions = [r for r in REGION_ORDER if r in all_df_ref['region'].unique()] if is_africa else []

    cat_data = []

    for (yi, yf), info in matrices.items():
        mat, dur = info["matrix"], info["duration"]
        if dur <= 0: continue

        row_sums = mat.sum(axis=1)  # Initial area per class
        col_sums = mat.sum(axis=0)  # Final area per class
        diag = pd.Series(np.diag(mat.values), index=mat.index)

        gross_gain = col_sums - diag
        gross_loss = row_sums - diag
        total_change = mat.values.sum() - diag.sum()
        S_t = (total_change / mat.values.sum()) / dur * 100 if mat.values.sum() > 0 else 0

        for c in lc_ids:
            G_int = (gross_gain[c] / dur) / col_sums[c] * 100 if col_sums[c] > 0 else 0
            L_int = (gross_loss[c] / dur) / row_sums[c] * 100 if row_sums[c] > 0 else 0

            row_data = {
                'interval': f"{yi}-{yf}",
                'year_initial': yi,
                'year_final': yf,
                'class_id': c,
                'class_name': LC_CLASSES[c],
                'gain_area_km2': round(gross_gain[c], 2),
                'loss_area_km2': round(gross_loss[c], 2),
                'gain_annual_km2_yr': round(gross_gain[c] / dur, 2),
                'loss_annual_km2_yr': round(gross_loss[c] / dur, 2),
                'gain_intensity_pct_yr': round(G_int, 4),
                'loss_intensity_pct_yr': round(L_int, 4),
                'uniform_intensity': round(S_t, 4)
            }

            # Add regional breakdown for AFRICA
            if is_africa:
                for reg in regions:
                    reg_mat = get_region_matrix(all_df_ref, reg, yi, yf)

                    # Regional gain for class c
                    reg_gain = reg_mat[c].sum() - reg_mat.at[c, c]
                    # Regional loss for class c
                    reg_loss = reg_mat.loc[c].sum() - reg_mat.at[c, c]

                    g_pct = (reg_gain / gross_gain[c] * 100) if gross_gain[c] > 0 else 0
                    l_pct = (reg_loss / gross_loss[c] * 100) if gross_loss[c] > 0 else 0

                    row_data[f'{reg}_gain_km2'] = round(reg_gain, 2)
                    row_data[f'{reg}_gain_pct'] = round(g_pct, 2)
                    row_data[f'{reg}_loss_km2'] = round(reg_loss, 2)
                    row_data[f'{reg}_loss_pct'] = round(l_pct, 2)

            cat_data.append(row_data)

    df_cat = pd.DataFrame(cat_data).sort_values(['year_initial', 'class_id'])
    if df_cat.empty: return
    save_table(df_cat, output_folder, f"L2_Category_{site_name}")

    # Plotting - MODIFIED for merged panels at zero
    intervals = df_cat['interval'].unique()
    n_int = len(intervals)

    # INCREASED figure height for better visibility
    fig_h = n_int * 4.5 + 4.0

    # Create subplots: Area (left) | Intensity (right) | Regional (if AFRICA)
    # Using 3 columns for merged view: Area, Intensity, Regional(optional)
    if is_africa:
        n_cols = 3
        w_ratios = [1, 1, 0.8]
        fig_w = 18
    else:
        n_cols = 2
        w_ratios = [1, 1]
        fig_w = 14

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Create GridSpec with different spacing:
    # - No space between Area and Intensity (wspace=0 for first two columns)
    # - Space before Regional panel
    from matplotlib.gridspec import GridSpec

    if is_africa:
        # Use nested gridspec for different spacing
        gs_main = GridSpec(n_int, 2, figure=fig, wspace=0.15, hspace=0.45,
                           width_ratios=[2, 0.8])  # [merged Area+Int, Regional]

        axes = []
        for i in range(n_int):
            # Create sub-gridspec for Area and Intensity (no space between them)
            gs_sub = gs_main[i, 0].subgridspec(1, 2, wspace=0, width_ratios=[1, 1])
            ax_area = fig.add_subplot(gs_sub[0, 0])
            ax_int = fig.add_subplot(gs_sub[0, 1])
            ax_reg = fig.add_subplot(gs_main[i, 1])
            axes.append({'area': ax_area, 'int': ax_int, 'reg': ax_reg})
    else:
        gs_main = GridSpec(n_int, 2, figure=fig, wspace=0, hspace=0.45,
                           width_ratios=[1, 1])
        axes = []
        for i in range(n_int):
            ax_area = fig.add_subplot(gs_main[i, 0])
            ax_int = fig.add_subplot(gs_main[i, 1])
            axes.append({'area': ax_area, 'int': ax_int})

    y = np.arange(len(classes))
    bh = 0.40  # Bar height
    # Gain bar at y + bh/2, Loss bar at y - bh/2
    # This makes them touch at y (no overlap, no gap between gain/loss)
    # Gap between categories = 1.0 - 2*bh = 1.0 - 0.8 = 0.2

    for i, interval in enumerate(intervals):
        data = df_cat[df_cat['interval'] == interval].set_index('class_name').reindex(classes)

        # Colors for panels 1 & 2: either LC colors or grayscale
        if grayscale_bars:
            cols = ['#696969'] * len(classes)  # Dark gray for all
            edge_cols = ['#404040'] * len(classes)  # Darker gray edges
        else:
            cols = [LC_COLORS[c] for c in classes]
            edge_cols = cols

        yi_val, yf_val = int(data['year_initial'].iloc[0]), int(data['year_final'].iloc[0])
        U = data['uniform_intensity'].iloc[0]

        ax_area = axes[i]['area']
        ax_int = axes[i]['int']
        if is_africa:
            ax_reg = axes[i]['reg']

        # =====================================================
        # PANEL 1 (LEFT): AREA (km²/yr) - bars go LEFT from zero
        # =====================================================
        ax_area.barh(y + bh/2, data['gain_annual_km2_yr'], height=bh,
                     color=cols, edgecolor='black', lw=0.5)
        ax_area.barh(y - bh/2, data['loss_annual_km2_yr'], height=bh,
                     color='white', edgecolor=edge_cols, hatch='////', lw=0.5)

        # Invert x-axis so bars go left, but zero is on the RIGHT edge
        ax_area.invert_xaxis()

        # Set x-axis to start from 0 on the right
        max_area = max(data['gain_annual_km2_yr'].max(), data['loss_annual_km2_yr'].max()) * 1.15
        ax_area.set_xlim(max_area, 0)  # Inverted: max on left, 0 on right

        # =====================================================
        # PANEL 2 (RIGHT): INTENSITY (%/yr) - bars go RIGHT from zero
        # =====================================================
        ax_int.barh(y + bh/2, data['gain_intensity_pct_yr'], height=bh,
                    color=cols, edgecolor='black', lw=0.5)
        ax_int.barh(y - bh/2, data['loss_intensity_pct_yr'], height=bh,
                    color='white', edgecolor=edge_cols, hatch='////', lw=0.5)

        # Uniform intensity line
        ax_int.axvline(U, color='red', linestyle='--', lw=1.5, zorder=10)

        # Set x-axis to start from 0 on the left
        max_int = max(data['gain_intensity_pct_yr'].max(), data['loss_intensity_pct_yr'].max(), U) * 1.2
        ax_int.set_xlim(0, max_int)

        # =====================================================
        # SHARED ZERO AXIS STYLING
        # =====================================================

        # Set consistent y-limits for all panels
        y_min, y_max = -0.5, len(classes) - 0.5
        ax_area.set_ylim(y_min, y_max)
        ax_int.set_ylim(y_min, y_max)
        if is_africa:
            ax_reg.set_ylim(y_min, y_max)

        # Add U value annotation at top of the red line (after y_max is defined)
        ax_int.text(U, y_max + 0.3, f'U={U:.2f}', ha='center', va='bottom',
                    fontsize=10, color='red', fontweight='bold', clip_on=False)

        # Remove inner spines where panels meet
        ax_area.spines['right'].set_linewidth(2)
        ax_area.spines['right'].set_color('black')
        ax_int.spines['left'].set_visible(False)

        # Y-axis: only on the left panel (Area)
        ax_area.set_yticks(y)
        ax_area.set_yticklabels(classes, fontsize=12, fontweight='bold')
        ax_area.set_ylabel(interval, fontweight='bold', fontsize=13)
        ax_area.tick_params(axis='y', which='both', left=True, labelleft=True)

        # Hide y-axis ticks on other panels
        ax_int.set_yticks([])
        ax_int.tick_params(axis='y', which='both', left=False, labelleft=False)
        if is_africa:
            ax_reg.set_yticks([])
            ax_reg.tick_params(axis='y', which='both', left=False, labelleft=False)

        # X-axis formatting
        ax_area.xaxis.set_major_formatter(ticker.FuncFormatter(human_format_abs))
        ax_area.tick_params(axis='x', labelsize=11)
        ax_int.tick_params(axis='x', labelsize=11)

        # Remove zero from both axes and add single "0" at the meeting point
        area_ticks = ax_area.get_xticks()
        int_ticks = ax_int.get_xticks()

        # Filter out zero ticks
        ax_area.set_xticks([t for t in area_ticks if t > 0.001])
        ax_int.set_xticks([t for t in int_ticks if t > 0.001])

        # Add single "0" label at the boundary between panels
        # Position it at the right edge of ax_area (which is where zero is)
        ax_area.text(0, -0.6, '0', ha='center', va='top', fontsize=11, fontweight='bold',
                     transform=ax_area.get_xaxis_transform(), clip_on=False)

        # Remove top spines
        ax_area.spines['top'].set_visible(False)
        ax_int.spines['top'].set_visible(False)

        # Grid lines
        ax_area.grid(axis='x', linestyle='--', alpha=0.3)
        ax_int.grid(axis='x', linestyle='--', alpha=0.3)

        # =====================================================
        # PANEL 3: REGIONAL CONTRIBUTION (%) [AFRICA only]
        # =====================================================
        if is_africa:
            ax_int.spines['right'].set_linewidth(1.5)
            ax_reg.spines['left'].set_visible(False)
            ax_reg.spines['top'].set_visible(False)

            # For each class, get regional % of gain and loss
            gain_contribs = {r: [] for r in regions}
            loss_contribs = {r: [] for r in regions}

            for c_name in classes:
                for reg in regions:
                    g_pct = data.loc[c_name, f'{reg}_gain_pct'] if f'{reg}_gain_pct' in data.columns else 0
                    l_pct = data.loc[c_name, f'{reg}_loss_pct'] if f'{reg}_loss_pct' in data.columns else 0
                    gain_contribs[reg].append(g_pct)
                    loss_contribs[reg].append(l_pct)

            # Plot stacked bars
            left_g = np.zeros(len(classes))
            left_l = np.zeros(len(classes))

            for reg in regions:
                ax_reg.barh(y + bh/2, gain_contribs[reg], height=bh, left=left_g,
                           color=REGION_COLORS[reg], edgecolor='white', lw=0.3)
                left_g += np.array(gain_contribs[reg])

                ax_reg.barh(y - bh/2, loss_contribs[reg], height=bh, left=left_l,
                           color=REGION_COLORS[reg], edgecolor='white', lw=0.3,
                           hatch='////', alpha=0.7)
                left_l += np.array(loss_contribs[reg])

            ax_reg.set_xlim(0, 100)
            ax_reg.tick_params(axis='x', labelsize=11)

    # Column titles (only on first row)
    axes[0]['area'].set_title("← Area (km²/yr)", fontsize=12, fontweight='bold', loc='left')
    axes[0]['int'].set_title("Intensity (%/yr) →", fontsize=12, fontweight='bold', loc='right')
    if is_africa:
        axes[0]['reg'].set_title("Regional %", fontsize=12, fontweight='bold')

    fig.suptitle(f"{site_name} - Category Level Intensity Analysis", fontsize=15, fontweight='bold', y=0.99)

    # Legend
    legend_items = [
        patches.Patch(facecolor='gray', edgecolor='black', label='Gain'),
        patches.Patch(facecolor='white', edgecolor='gray', hatch='////', label='Loss'),
        mlines.Line2D([0], [0], color='red', linestyle='--', label='Uniform')
    ]
    if is_africa:
        for reg in regions:
            legend_items.append(patches.Patch(facecolor=REGION_COLORS[reg], label=reg))

    fig.legend(handles=legend_items, loc='lower center', bbox_to_anchor=(0.5, 0.01),
              ncol=len(legend_items), frameon=False, fontsize=11)

    plt.subplots_adjust(left=0.08, right=0.98, bottom=0.07, top=0.95)
    suffix = "_grayscale" if grayscale_bars else ""
    plt.savefig(os.path.join(output_folder, f"L2_Category_{site_name}{suffix}.png"), bbox_inches='tight')
    plt.close()

# -------------------------------------------------------------------
# 5. TRANSITION LEVEL (Standard IA - unchanged structure)
# -------------------------------------------------------------------

def analyze_transition_level(matrices, site_name, all_df_ref=None):
    print(f"   [Level 3] Transition Analysis: {site_name}")

    lc_ids = sorted(LC_CLASSES.keys())
    n_trans = len(lc_ids) * (len(lc_ids) - 1)
    results = []
    intervals_list = []

    is_africa = (site_name == "AFRICA") and (all_df_ref is not None)
    regions = [r for r in REGION_ORDER if r in all_df_ref['region'].unique()] if is_africa else []

    for (yi, yf), info in matrices.items():
        interval = f"{yi}-{yf}"
        if interval not in intervals_list: intervals_list.append(interval)

        mat, dur = info["matrix"], info["duration"]
        if dur <= 0: continue

        initial = mat.sum(axis=1)
        diag = np.diag(mat.values)
        total_change = mat.values.sum() - diag.sum()
        uniform_qty = total_change / n_trans if n_trans > 0 else 0

        for n in lc_ids:  # Gaining
            gross_gain_n = mat[n].sum() - mat.at[n, n]
            other_area = initial.sum() - initial[n]
            W_tn = (gross_gain_n / dur) / other_area * 100 if other_area > 0 else 0

            for m in lc_ids:  # Losing
                if m == n: continue
                trans_area = mat.at[m, n]
                R_tin = (trans_area / dur) / initial[m] * 100 if initial[m] > 0 else 0

                row_data = {
                    'interval': interval,
                    'year_initial': yi,
                    'year_final': yf,
                    'from_class': LC_CLASSES[m],
                    'to_class': LC_CLASSES[n],
                    'from_id': m,
                    'to_id': n,
                    'transition': f"{LC_CLASSES[m]}→{LC_CLASSES[n]}",
                    'area_km2': round(trans_area, 2),
                    'R_tin': round(R_tin, 4),
                    'W_tn': round(W_tn, 4),
                    'gain_diff': round(R_tin - W_tn, 4),
                    'is_targeted': R_tin > W_tn,
                    'uniform_qty_threshold': round(uniform_qty, 2)
                }

                # Add regional breakdown for AFRICA
                if is_africa and trans_area > 0:
                    for reg in regions:
                        reg_mat = get_region_matrix(all_df_ref, reg, yi, yf)
                        reg_trans = reg_mat.at[m, n]
                        reg_pct = (reg_trans / trans_area * 100) if trans_area > 0 else 0

                        row_data[f'{reg}_area_km2'] = round(reg_trans, 2)
                        row_data[f'{reg}_contrib_pct'] = round(reg_pct, 2)

                results.append(row_data)

    df_trans = pd.DataFrame(results)
    if df_trans.empty: return
    save_table(df_trans, output_folder, f"L3_Transition_{site_name}")

    # Plot stationarity matrix (with regional decomposition for AFRICA)
    plot_stationarity_matrix(df_trans, intervals_list, site_name, all_df_ref)

    # Plot heatmap grids
    plot_transition_heatmaps(matrices, site_name)

def plot_stationarity_matrix(df, intervals, site_name, all_df_ref=None):
    classes = [LC_CLASSES[i] for i in sorted(LC_CLASSES.keys())]
    N = len(classes)
    n_int = len(intervals)
    is_africa = (site_name == "AFRICA") and (all_df_ref is not None)

    fig, ax = plt.subplots(figsize=(14, 14))
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)

    # Pre-calculate regional contributions for AFRICA
    if is_africa:
        regions = [r for r in REGION_ORDER if r in all_df_ref['region'].unique()]

        # Build lookup: (interval, from_id, to_id) -> {region: area}
        reg_contrib_cache = {}
        for interval in intervals:
            yi, yf = map(int, interval.split('-'))
            for reg in regions:
                reg_mat = get_region_matrix(all_df_ref, reg, yi, yf)
                for m in sorted(LC_CLASSES.keys()):
                    for n in sorted(LC_CLASSES.keys()):
                        if m == n: continue
                        key = (interval, m, n)
                        if key not in reg_contrib_cache:
                            reg_contrib_cache[key] = {}
                        reg_contrib_cache[key][reg] = reg_mat.at[m, n]

    for r_idx, from_cls in enumerate(classes):
        for c_idx, to_cls in enumerate(classes):
            y_base = N - 1 - r_idx
            x_base = c_idx

            # Diagonal = Persistence
            if from_cls == to_cls:
                ax.add_patch(patches.Rectangle((x_base, y_base), 1, 1,
                            facecolor='#dddddd', hatch='///'))
                continue

            strip_w = 1.0 / n_int
            cell = df[(df['from_class'] == from_cls) & (df['to_class'] == to_cls)]

            # Build info map for this cell
            info_map = {}
            for _, r in cell.iterrows():
                info_map[r['interval']] = {
                    'targeted': r['is_targeted'],
                    'area': r['area_km2'],
                    'from_id': r['from_id'],
                    'to_id': r['to_id']
                }

            for k, intv in enumerate(intervals):
                info = info_map.get(intv, {'targeted': False, 'area': 0})
                targeted = info['targeted']
                rect_x = x_base + k * strip_w

                # For AFRICA: show regional decomposition for targeted transitions
                if is_africa and targeted and info['area'] > 0:
                    # Get regional contributions
                    key = (intv, info['from_id'], info['to_id'])
                    reg_areas = reg_contrib_cache.get(key, {})
                    total_area = sum(reg_areas.values())

                    if total_area > 0:
                        # Draw stacked rectangles from bottom to top
                        y_offset = 0
                        for reg in regions:
                            reg_area = reg_areas.get(reg, 0)
                            height_frac = reg_area / total_area
                            if height_frac > 0.001:  # Skip tiny slivers
                                ax.add_patch(patches.Rectangle(
                                    (rect_x, y_base + y_offset),
                                    strip_w, height_frac,
                                    facecolor=REGION_COLORS[reg],
                                    edgecolor='none', lw=0
                                ))
                            y_offset += height_frac
                    else:
                        # Fallback to gray
                        ax.add_patch(patches.Rectangle((rect_x, y_base), strip_w, 1,
                                    facecolor='#f0f0f0', edgecolor='none', lw=0))

                # Non-AFRICA: solid colors based on targeted only
                elif not is_africa:
                    if targeted:
                        color = "#d62728"  # Red: Targeted
                    else:
                        color = "#f0f0f0"  # Gray: Avoided

                    ax.add_patch(patches.Rectangle((rect_x, y_base), strip_w, 1,
                                facecolor=color, edgecolor='none', lw=0))

                # AFRICA but not targeted: gray
                else:
                    ax.add_patch(patches.Rectangle((rect_x, y_base), strip_w, 1,
                                facecolor='#f0f0f0', edgecolor='none', lw=0))

                # Draw black separator line after each strip (except last)
                if k < n_int - 1:
                    line_x = rect_x + strip_w
                    ax.plot([line_x, line_x], [y_base, y_base + 1],
                           color='black', lw=0.3, zorder=10)

            # Cell border - THICKER to distinguish from interval separators
            ax.add_patch(patches.Rectangle((x_base, y_base), 1, 1,
                        fill=False, edgecolor='black', lw=2.5))

    # Grid lines - THICKER main category separators
    for i in range(N + 1):
        ax.axhline(i, color='black', lw=2.5)
        ax.axvline(i, color='black', lw=2.5)

    ax.set_xticks(np.arange(N) + 0.5)
    ax.set_yticks(np.arange(N) + 0.5)
    ax.set_xticklabels(classes, rotation=90, fontsize=14, fontweight='bold')
    ax.set_yticklabels(classes[::-1], fontsize=14, fontweight='bold')
    ax.set_xlabel("To Category (Gaining)", fontsize=16, fontweight='bold')
    ax.set_ylabel("From Category (Losing)", fontsize=16, fontweight='bold')

    if is_africa:
        ax.set_title(f"{site_name}: Transition-Level Stationarity with Regional Attribution\n(Strips: {intervals[0]} → {intervals[-1]})",
                    fontsize=16, fontweight='bold')
        # Legend with regions - positioned lower and larger
        legend = [
            patches.Patch(facecolor='#f0f0f0', edgecolor='black', label='Avoided/Random'),
            patches.Patch(facecolor='#dddddd', hatch='///', edgecolor='black', label='Persistence')
        ]
        for reg in regions:
            legend.append(patches.Patch(facecolor=REGION_COLORS[reg], label=reg))
        ax.legend(handles=legend, loc='upper center', bbox_to_anchor=(0.5, -0.10),
                 ncol=4, frameon=False, fontsize=13)
    else:
        ax.set_title(f"{site_name}: Transition-Level Stationarity\n(Strips: {intervals[0]} → {intervals[-1]})",
                    fontsize=16, fontweight='bold')
        legend = [
            patches.Patch(facecolor='#d62728', edgecolor='black', label='Targeted'),
            patches.Patch(facecolor='#f0f0f0', edgecolor='black', label='Avoided/Random'),
            patches.Patch(facecolor='#dddddd', hatch='///', edgecolor='black', label='Persistence')
        ]
        ax.legend(handles=legend, loc='upper center', bbox_to_anchor=(0.5, -0.10),
                 ncol=3, frameon=False, fontsize=13)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(os.path.join(output_folder, f"L3_Stationarity_{site_name}.png"))
    plt.close()

def plot_transition_heatmaps(matrices, site_name):
    intervals = list(matrices.keys())
    n = len(intervals)
    if n == 0: return

    cols = 2
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(11*cols, 9*rows), constrained_layout=True)
    axes = np.array(axes).flatten() if n > 1 else [axes]

    labels = [LC_CLASSES[i] for i in sorted(LC_CLASSES.keys())]

    for i, (yi, yf) in enumerate(intervals):
        ax = axes[i]
        mat = matrices[(yi, yf)]["matrix"]
        pct = mat.div(mat.sum(axis=1), axis=0).fillna(0) * 100

        # Annotation: percentage + raw area
        annot = pd.DataFrame(index=mat.index, columns=mat.columns)
        for r in mat.index:
            for c in mat.columns:
                annot.at[r, c] = f"{pct.at[r,c]:.1f}%\n({mat.at[r,c]:.0f})"

        sns.heatmap(pct, ax=ax, annot=annot, fmt="", cmap="OrRd", vmin=0, vmax=100,
                   cbar=False, xticklabels=labels, yticklabels=labels,
                   annot_kws={"size": 10, "weight": "bold"}, linewidths=1, linecolor='gray')

        ax.set_title(f"{yi}-{yf}", fontsize=16, fontweight='bold')
        ax.set_ylabel("From", fontsize=13, fontweight='bold')
        ax.set_xlabel("To", fontsize=13, fontweight='bold')
        ax.tick_params(axis='both', labelsize=12)
        for tick in ax.get_xticklabels(): tick.set_fontweight('bold')
        for tick in ax.get_yticklabels(): tick.set_fontweight('bold')

    for j in range(i + 1, len(axes)): axes[j].axis('off')

    norm = plt.Normalize(0, 100)
    sm = plt.cm.ScalarMappable(cmap="OrRd", norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='horizontal', fraction=0.04, pad=0.02, aspect=50)
    cbar.set_label('Row Fraction (%)', fontsize=14, fontweight='bold')

    fig.suptitle(f"{site_name}: Transition Matrices", fontsize=18, fontweight='bold')
    plt.savefig(os.path.join(output_folder, f"L3_Heatmaps_{site_name}.png"))
    plt.close()

# -------------------------------------------------------------------
# 7. CROSS-REGIONAL COMPARISON FIGURE
# -------------------------------------------------------------------

def plot_cross_regional_comparison(all_df):
    """
    Single figure comparing all regions' intensities side-by-side.
    This shows RATE differences, complementing the AFRICA decomposition (SIZE).
    """
    print("\n   [Cross-Regional] Building comparison figure...")

    regions = [r for r in REGION_ORDER if r in all_df['region'].unique()]

    # Collect metrics for all regions
    all_region_metrics = {}
    for reg in regions:
        df_site = all_df[all_df['region'] == reg]
        matrices = build_matrices(df_site)
        metrics, U = compute_interval_metrics(matrices)
        metrics['region'] = reg
        metrics['U'] = U
        all_region_metrics[reg] = metrics

    # Get Africa metrics for reference
    df_africa = all_df[all_df['region'] == 'AFRICA']
    africa_matrices = build_matrices(df_africa)
    africa_metrics, africa_U = compute_interval_metrics(africa_matrices)

    intervals = africa_metrics['interval'].tolist()

    # Create figure with 2 panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'width_ratios': [1.2, 1]})

    # --- PANEL 1: Grouped Bar Chart - Intensity by Region and Interval ---
    x = np.arange(len(intervals))
    n_regions = len(regions)
    width = 0.8 / n_regions

    for i, reg in enumerate(regions):
        reg_metrics = all_region_metrics[reg].set_index('interval')
        intensities = [reg_metrics.loc[intv, 'intensity'] if intv in reg_metrics.index else 0
                      for intv in intervals]
        offset = (i - n_regions/2 + 0.5) * width
        ax1.bar(x + offset, intensities, width=width, color=REGION_COLORS[reg],
               edgecolor='black', linewidth=0.5, label=reg)

    # Africa uniform line
    ax1.axhline(africa_U, color='black', linestyle='--', linewidth=2, label=f'Africa Uniform: {africa_U:.2f}%')

    ax1.set_xticks(x)
    ax1.set_xticklabels(intervals, rotation=45, ha='right', fontsize=11, fontweight='bold')
    ax1.set_ylabel("Annual Change Intensity (%/yr)", fontsize=12, fontweight='bold')
    ax1.set_xlabel("Time Interval", fontsize=12, fontweight='bold')
    ax1.set_title("Regional Intensity Comparison by Interval", fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(axis='y', linestyle='--', alpha=0.5)

    # --- PANEL 2: Summary Statistics ---
    summary_data = []
    for reg in regions:
        reg_metrics = all_region_metrics[reg]
        mean_int = reg_metrics['intensity'].mean()
        max_int = reg_metrics['intensity'].max()
        max_interval = reg_metrics.loc[reg_metrics['intensity'].idxmax(), 'interval']
        total_change = reg_metrics['change_area'].sum()

        # Deviation from Africa uniform
        deviation = ((mean_int - africa_U) / africa_U) * 100

        summary_data.append({
            'Region': reg,
            'Mean Intensity': mean_int,
            'Peak Intensity': max_int,
            'Peak Interval': max_interval,
            'Total Change (km²)': total_change,
            'vs Africa (%)': deviation
        })

    df_summary = pd.DataFrame(summary_data)

    # Save detailed summary with all intervals
    detailed_rows = []
    for reg in regions:
        reg_metrics = all_region_metrics[reg]
        for _, row in reg_metrics.iterrows():
            detailed_rows.append({
                'region': reg,
                'interval': row['interval'],
                'year_initial': row['year_initial'],
                'year_final': row['year_final'],
                'duration_yrs': row['duration'],
                'total_area_km2': round(row['total_area'], 2),
                'change_area_km2': round(row['change_area'], 2),
                'intensity_pct_yr': round(row['intensity'], 4),
                'change_pct_of_map': round(row['change_pct'], 4)
            })

    df_detailed = pd.DataFrame(detailed_rows)
    save_table(df_detailed, output_folder, "Cross_Regional_AllIntervals")

    # Horizontal bar chart of mean intensities
    y_pos = np.arange(len(regions))
    colors = [REGION_COLORS[r] for r in regions]

    bars = ax2.barh(y_pos, df_summary['Mean Intensity'], color=colors, edgecolor='black', height=0.6)
    ax2.axvline(africa_U, color='black', linestyle='--', linewidth=2, label=f'Africa Uniform')

    # Add value labels
    for i, bar in enumerate(bars):
        mean_int = df_summary.iloc[i]['Mean Intensity']
        vs_africa = df_summary.iloc[i]['vs Africa (%)']
        sign = '+' if vs_africa > 0 else ''
        ax2.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height()/2,
                f'{mean_int:.2f}% ({sign}{vs_africa:.0f}%)', va='center', fontsize=10, fontweight='bold')

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(regions, fontsize=12, fontweight='bold')
    ax2.set_xlabel("Mean Annual Intensity (%/yr)", fontsize=12, fontweight='bold')
    ax2.set_title("Regional Mean Intensity\n(values show deviation from Africa)", fontsize=14, fontweight='bold')
    ax2.grid(axis='x', linestyle='--', alpha=0.5)

    # Extend x-axis for labels
    xlim = ax2.get_xlim()
    ax2.set_xlim(xlim[0], xlim[1] * 1.4)

    fig.suptitle("Cross-Regional Intensity Comparison", fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)

    plt.savefig(os.path.join(output_folder, "Cross_Regional_Comparison.png"), bbox_inches='tight')
    plt.close()

    # Save summary table
    save_table(df_summary, output_folder, "Cross_Regional_Summary")
    print("   Created: Cross_Regional_Comparison.png")

# -------------------------------------------------------------------
# 6. MAIN
# -------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("INTENSITY ANALYSIS WITH REGIONAL DECOMPOSITION")
    print("=" * 60 + "\n")

    # Process ONLY AFRICA (main figures with decomposition)
    print("\nProcessing: AFRICA (with regional decomposition)")
    df_africa = all_df[all_df['region'] == 'AFRICA']
    matrices = build_matrices(df_africa)

    # Level 1: Interval
    # Use compute_interval_metrics for plotting
    df_int, U = compute_interval_metrics(matrices)
    # Use compute_interval_with_regional for CSV (includes 3rd panel data)
    df_int_with_regional, _ = compute_interval_with_regional(matrices, all_df)
    save_table(df_int_with_regional, output_folder, "L1_Interval_AFRICA")
    plot_interval_level(df_int, U, "AFRICA", all_df)

    # Level 2: Category (MODIFIED - merged panels)
    analyze_category_level(matrices, "AFRICA", all_df, grayscale_bars=False)  # Color version
    analyze_category_level(matrices, "AFRICA", all_df, grayscale_bars=True)   # Grayscale version

    # Level 3: Transition
    analyze_transition_level(matrices, "AFRICA", all_df)

    # Cross-Regional Comparison
    plot_cross_regional_comparison(all_df)

    print(f"\n{'=' * 60}")
    print(f"Complete! Output: {output_folder}")
    print("=" * 60)
    print("\nFigures for paper:")
    print("  1. L1_Interval_AFRICA.png          - Interval level + regional contribution")
    print("  2. L2_Category_AFRICA.png          - Category level (color) + regional contribution")
    print("  3. L2_Category_AFRICA_grayscale.png - Category level (grayscale panels 1&2)")
    print("  4. L3_Stationarity_AFRICA.png      - Transition targeting + regional attribution")
    print("  5. L3_Heatmaps_AFRICA.png          - Raw transition matrices")
    print("  6. Cross_Regional_Comparison.png   - Regional intensity comparison")
