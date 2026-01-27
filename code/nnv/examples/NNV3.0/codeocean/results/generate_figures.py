#!/usr/bin/env python3
"""Generate publication-quality figures and LaTeX tables from CodeOcean results.

Outputs (saved to figures/):
  - fig_weight_pert_mnist_mlp.pdf/png  (ModelStar bar chart)
  - fig_gnnv_results.pdf/png           (GNNV sensitivity line chart)
  - fig_fairness_individual.pdf/png    (FairNNV stacked area plot)
  - tab_fairness_results.tex           (FairNNV LaTeX table)
  - tab_videostar_results.tex          (VideoStar LaTeX table)
  - fig_probver_results.pdf/png        (ProbVer horizontal bar chart)
  - tab_probver_results.tex            (ProbVer LaTeX table)
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import scipy.io

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
RESULTS_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(RESULTS_DIR, 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

DPI = 300

# Publication style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

# Colorblind-friendly palette
BLUE   = (0.0, 0.45, 0.7)
ORANGE = (0.9, 0.6,  0.0)
TEAL   = (0.0, 0.6,  0.5)
GREEN  = (46/255, 204/255, 113/255)
RED    = (231/255, 76/255, 60/255)


# ---------------------------------------------------------------------------
# 1. ModelStar: Weight Perturbation Verification
# ---------------------------------------------------------------------------
def generate_modelstar_figure():
    csv_path = os.path.join(RESULTS_DIR, 'modelstar_results.csv')
    df = pd.read_csv(csv_path)

    layers = df['Layer'].unique()
    n_layers = len(layers)
    fig, axes = plt.subplots(2, n_layers, figsize=(3.2 * n_layers, 4.8),
                             gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.35})
    if n_layers == 1:
        axes = axes.reshape(2, 1)

    layer_labels = {'fc_6': 'Layer fc_6', 'fc_5': 'Layer fc_5', 'fc_4': 'Layer fc_4'}

    for col, layer in enumerate(layers):
        layer_df = df[df['Layer'] == layer]
        fracs = layer_df['PerturbationFrac'].values
        verified = layer_df['VerifiedPercent'].values
        avg_time = layer_df['AvgTimePerImage'].values

        # Top row: verified %
        ax_top = axes[0, col]
        bars = ax_top.bar(range(len(fracs)), verified, color=BLUE, edgecolor='white', width=0.6)
        for bar, v in zip(bars, verified):
            if v > 5:
                ax_top.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 3,
                            f'{v:.0f}', ha='center', va='top', fontsize=8, color='white', fontweight='bold')
            else:
                ax_top.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                            f'{v:.0f}', ha='center', va='bottom', fontsize=8, color='black')
        ax_top.set_xticks(range(len(fracs)))
        ax_top.set_xticklabels([f'{f*100:.1f}%' for f in fracs])
        ax_top.set_title(layer_labels.get(layer, layer))
        ax_top.set_ylim(0, 105)
        if col == 0:
            ax_top.set_ylabel('Verified (%)')

        # Bottom row: avg time per image
        ax_bot = axes[1, col]
        bars_t = ax_bot.bar(range(len(fracs)), avg_time, color=ORANGE, edgecolor='white', width=0.6)
        for bar, t in zip(bars_t, avg_time):
            ax_bot.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(avg_time)*0.02,
                        f'{t:.2f}', ha='center', va='bottom', fontsize=7)
        ax_bot.set_xticks(range(len(fracs)))
        ax_bot.set_xticklabels([f'{f*100:.1f}%' for f in fracs])
        ax_bot.set_xlabel('Perturbation Fraction')
        if col == 0:
            ax_bot.set_ylabel('Avg. Time / Image (s)')

    # Share y-axis on top row only (verified % is always 0-100)
    for col in range(1, n_layers):
        axes[0, col].sharey(axes[0, 0])
        axes[0, col].set_yticklabels([])

    fig.suptitle('ModelStar Weight Perturbation Verification (MNIST MLP)', fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    save_figure(fig, 'fig_weight_pert_mnist_mlp')
    print('  Generated: fig_weight_pert_mnist_mlp')


# ---------------------------------------------------------------------------
# 2. GNNV: GNN Verification Sensitivity
# ---------------------------------------------------------------------------
def generate_gnnv_figure():
    mat_path = os.path.join(RESULTS_DIR, 'gnn_results.mat')
    mat = scipy.io.loadmat(mat_path, squeeze_me=False)

    results = mat['results'][0, 0]
    config = results['config'][0, 0]

    # Extract config
    models = [str(m).strip().strip("[]'\"") for m in config[0].flatten()]
    epsilons = config[1].flatten().astype(float)
    num_scenarios = int(config[5].flatten()[0])

    data = results['data']  # shape: (n_models, n_eps, n_scenarios) of structs

    # Total voltage nodes (nodes with status != -1)
    first_entry = data[0, 0, 0][0, 0]
    verif_per_node = first_entry['verif_per_node'].flatten()
    total_voltage_nodes = np.sum(verif_per_node != -1)

    # Compute verified % for each (model, epsilon)
    pct = np.zeros((len(models), len(epsilons)))
    for m in range(len(models)):
        for e in range(len(epsilons)):
            total_verified = 0
            for s in range(num_scenarios):
                entry = data[m, e, s][0, 0]
                total_verified += int(entry['verified'].flatten()[0])
            pct[m, e] = (total_verified / (total_voltage_nodes * num_scenarios)) * 100

    # Plot
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    markers = ['o', 's', 'D']
    colors = [BLUE, ORANGE, TEAL]

    for m, (model, marker, color) in enumerate(zip(models, markers, colors)):
        ax.plot(epsilons, pct[m], marker=marker, color=color, linewidth=2.0,
                markersize=7, markerfacecolor=color, label=model)

    ax.set_xlabel(r'Perturbation $\epsilon$')
    ax.set_ylabel('Verified (%)')
    ax.set_title('GNN Verification vs. Perturbation (IEEE 24-Bus)')
    ax.set_xlim(0, max(epsilons) * 1.15)
    ax.set_ylim(0, 100)
    ax.legend(loc='lower left')

    plt.tight_layout()
    save_figure(fig, 'fig_gnnv_results')
    print('  Generated: fig_gnnv_results')


# ---------------------------------------------------------------------------
# 3. FairNNV: Individual Fairness Stacked Area Plot
# ---------------------------------------------------------------------------
def generate_fairness_figure():
    individual_csv = find_latest('fm26_individual_*.csv')
    df = pd.read_csv(individual_csv)

    models = df['Model'].unique()
    model_names = {'AC-1': 'Small Model (AC-1)', 'AC-3': 'Medium Model (AC-3)'}

    fig, axes = plt.subplots(1, len(models), figsize=(4.0 * len(models), 3.0), sharey=True)
    if len(models) == 1:
        axes = [axes]

    for ax, model in zip(axes, models):
        mdf = df[df['Model'] == model].sort_values('Epsilon')
        eps = mdf['Epsilon'].values
        fair = mdf['FairPercent'].values
        unfair = mdf['UnfairPercent'].values

        x = np.arange(len(eps))
        ax.fill_between(x, 0, unfair, color=RED, alpha=0.9, label='Unfair')
        ax.fill_between(x, unfair, unfair + fair, color=GREEN, alpha=0.9, label='Fair')
        ax.plot(x, unfair, 'w-', linewidth=1.5)

        ax.set_xticks(x)
        ax.set_xticklabels([f'{e:.2f}' for e in eps])
        ax.set_xlabel(r'Perturbation $\epsilon$')
        ax.set_title(model_names.get(model, model))
        ax.set_ylim(0, 100)
        ax.set_xlim(-0.3, len(eps) - 0.7)

    axes[0].set_ylabel('Percentage (%)')
    axes[-1].legend(loc='upper right')
    fig.suptitle('Individual Fairness Verification (Adult Census)', fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()

    save_figure(fig, 'fig_fairness_individual')
    print('  Generated: fig_fairness_individual')


# ---------------------------------------------------------------------------
# 4. FairNNV: LaTeX Table
# ---------------------------------------------------------------------------
def generate_fairness_table():
    counterfactual_csv = find_latest('fm26_counterfactual_*.csv')
    individual_csv = find_latest('fm26_individual_*.csv')
    timing_csv = find_latest('fm26_timing_*.csv')

    cf_df = pd.read_csv(counterfactual_csv)
    ind_df = pd.read_csv(individual_csv)
    tim_df = pd.read_csv(timing_csv)

    model_names = {'AC-1': 'Small (16-8)', 'AC-3': 'Medium (50)'}

    out_path = os.path.join(FIGURES_DIR, 'tab_fairness_results.tex')
    with open(out_path, 'w') as f:
        f.write('\\begin{table}[t]\n')
        f.write('\\centering\n')
        f.write('\\caption{Fairness Verification Results on Adult Census Dataset}\n')
        f.write('\\label{tab:fairness_results}\n')
        f.write('\\small\n')
        f.write('\\begin{tabular}{llrrr}\n')
        f.write('\\toprule\n')
        f.write('Model & $\\epsilon$ & VF (\\%) & Unfair (\\%) & Avg. Time (s) \\\\\n')
        f.write('\\midrule\n')

        # Counterfactual section
        f.write('\\multicolumn{5}{l}{\\textit{Counterfactual Fairness}} \\\\\n')
        for _, row in cf_df.iterrows():
            model = row['Model']
            display = model_names.get(model, model)
            # Find timing at eps=0
            t_row = tim_df[(tim_df['Model'] == model) & (tim_df['Epsilon'] == 0)]
            avg_time = t_row['AvgTimePerSample'].values[0] if len(t_row) else float('nan')
            f.write(f'{display} & 0.0 & {row["FairPercent"]:.0f} & {row["UnfairPercent"]:.0f} & {avg_time:.4f} \\\\\n')

        f.write('\\midrule\n')
        f.write('\\multicolumn{5}{l}{\\textit{Individual Fairness}} \\\\\n')

        for model in ind_df['Model'].unique():
            display = model_names.get(model, model)
            mdf = ind_df[ind_df['Model'] == model].sort_values('Epsilon')
            for i, (_, row) in enumerate(mdf.iterrows()):
                eps = row['Epsilon']
                t_row = tim_df[(tim_df['Model'] == model) & (np.isclose(tim_df['Epsilon'], eps))]
                avg_time = t_row['AvgTimePerSample'].values[0] if len(t_row) else float('nan')
                name = display if i == 0 else ''
                f.write(f'{name} & {eps:.2f} & {row["FairPercent"]:.0f} & {row["UnfairPercent"]:.0f} & {avg_time:.4f} \\\\\n')
            if model != ind_df['Model'].unique()[-1]:
                f.write('\\addlinespace\n')

        f.write('\\bottomrule\n')
        f.write('\\end{tabular}\n')
        f.write('\\end{table}\n')

    print('  Generated: tab_fairness_results.tex')


# ---------------------------------------------------------------------------
# 5. VideoStar: Parse RTF and generate LaTeX table
# ---------------------------------------------------------------------------
def generate_videostar_table():
    eps_values = ['1_255', '2_255', '3_255']
    eps_labels = ['1/255', '2/255', '3/255']
    rows = []

    for eps_val, eps_label in zip(eps_values, eps_labels):
        csv_path = os.path.join(RESULTS_DIR, f'eps={eps_val}.csv')
        data = parse_rtf_csv(csv_path)
        if data is not None:
            n_total = len(data)
            n_verified = sum(1 for r in data if r['Result'] == 1)
            n_unknown = sum(1 for r in data if r['Result'] == 2)
            avg_time = np.mean([r['Time'] for r in data])
            rows.append((eps_label, n_verified, n_unknown, n_total, avg_time))

    if not rows:
        print('  Skipped: VideoStar (no data)')
        return

    out_path = os.path.join(FIGURES_DIR, 'tab_videostar_results.tex')
    with open(out_path, 'w') as f:
        f.write('\\begin{table}[t]\n')
        f.write('\\centering\n')
        f.write('\\caption{VideoStar Verification Results (ZoomIn-4f)}\n')
        f.write('\\label{tab:videostar_results}\n')
        f.write('\\begin{tabular}{lcccc}\n')
        f.write('\\toprule\n')
        f.write('$\\epsilon$ & Verified & Unknown & Total & Avg. Time (s) \\\\\n')
        f.write('\\midrule\n')
        for eps_label, verified, unknown, total, avg_time in rows:
            f.write(f'${eps_label}$ & {verified} & {unknown} & {total} & {avg_time:.1f} \\\\\n')
        f.write('\\bottomrule\n')
        f.write('\\end{tabular}\n')
        f.write('\\end{table}\n')

    print('  Generated: tab_videostar_results.tex')


# ---------------------------------------------------------------------------
# 6. ProbVer: Probabilistic Verification Results
# ---------------------------------------------------------------------------
def generate_probver_figure():
    csv_path = os.path.join(RESULTS_DIR, 'results_summary.csv')
    df = pd.read_csv(csv_path)

    # Extract property ID from vnnlib filename
    df['property'] = df['vnnlib'].str.extract(r'prop_(\d+)')[0].astype(int)
    df['label'] = df['property'].apply(lambda p: f'Prop {p}')

    # Status colors
    status_colors = {'unsat': TEAL, 'sat': RED, 'unknown': ORANGE, 'error': (0.5, 0.5, 0.5)}

    fig, ax = plt.subplots(figsize=(5.0, 2.2))

    y_pos = np.arange(len(df))
    bar_colors = [status_colors.get(s, BLUE) for s in df['status']]

    bars = ax.barh(y_pos, df['time'], color=bar_colors, edgecolor='white', height=0.5)

    # Time labels on bars
    for bar, t, status in zip(bars, df['time'], df['status']):
        label = f'{t:.0f}s ({status.upper()})'
        ax.text(bar.get_width() - 5, bar.get_y() + bar.get_height()/2,
                label, ha='right', va='center', fontsize=8, color='white', fontweight='bold')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['label'])
    ax.set_xlabel('Verification Time (s)')
    ax.set_title('CP-Star Probabilistic Verification (TinyYOLO)')
    ax.invert_yaxis()
    ax.set_xlim(0, max(df['time']) * 1.05)

    plt.tight_layout()
    save_figure(fig, 'fig_probver_results')
    print('  Generated: fig_probver_results')


def generate_probver_table():
    csv_path = os.path.join(RESULTS_DIR, 'results_summary.csv')
    df = pd.read_csv(csv_path)

    df['property'] = df['vnnlib'].str.extract(r'prop_(\d+)')[0].astype(int)

    out_path = os.path.join(FIGURES_DIR, 'tab_probver_results.tex')
    with open(out_path, 'w') as f:
        f.write('\\begin{table}[t]\n')
        f.write('\\centering\n')
        f.write('\\caption{Probabilistic Verification Results (TinyYOLO, CP-Star)}\n')
        f.write('\\label{tab:probver_results}\n')
        f.write('\\begin{tabular}{lcrc}\n')
        f.write('\\toprule\n')
        f.write('Property & $\\epsilon$ & Time (s) & Result \\\\\n')
        f.write('\\midrule\n')
        for _, row in df.iterrows():
            prop_id = row['property']
            # Extract epsilon from filename
            eps_match = re.search(r'eps_(\d+)_(\d+)', row['vnnlib'])
            eps_str = f'{eps_match.group(1)}/{eps_match.group(2)}' if eps_match else '--'
            status = row['status'].upper()
            f.write(f'Prop {prop_id} & ${eps_str}$ & {row["time"]:.1f} & {status} \\\\\n')
        f.write('\\bottomrule\n')
        f.write('\\end{tabular}\n')
        f.write('\\end{table}\n')

    print('  Generated: tab_probver_results.tex')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def save_figure(fig, name):
    fig.savefig(os.path.join(FIGURES_DIR, f'{name}.pdf'), format='pdf')
    fig.savefig(os.path.join(FIGURES_DIR, f'{name}.png'), format='png')
    plt.close(fig)


def find_latest(pattern):
    files = glob.glob(os.path.join(RESULTS_DIR, pattern))
    if not files:
        raise FileNotFoundError(f'No files matching {pattern}')
    return max(files, key=os.path.getmtime)


def parse_rtf_csv(path):
    """Parse CSV data embedded in RTF files."""
    with open(path, 'r', errors='replace') as f:
        text = f.read()

    # Extract CSV lines from RTF
    lines = re.findall(r'(\d+,\d+,[\d.]+,\w+)', text)
    if not lines:
        return None

    data = []
    for line in lines:
        parts = line.split(',')
        data.append({
            'SampleNumber': int(parts[0]),
            'Result': int(parts[1]),
            'Time': float(parts[2]),
            'Method': parts[3],
        })
    return data


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print('Generating figures and tables from CodeOcean results...')
    print(f'  Results dir: {RESULTS_DIR}')
    print(f'  Output dir:  {FIGURES_DIR}')
    print()

    generate_modelstar_figure()
    generate_gnnv_figure()
    generate_fairness_figure()
    generate_fairness_table()
    generate_videostar_table()
    generate_probver_figure()
    generate_probver_table()

    print()
    print('Done. All outputs saved to figures/')
