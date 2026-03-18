#!/usr/bin/env python3
# ================================================================
#  branch_dos.py  — Cleaned version
#  Branch-resolved phonon DOS from QE matdyn .freq file
#  Gaussian smearing, vectorized (fast), works for any nbnd
# ================================================================

import re
import numpy as np
import matplotlib.pyplot as plt
import csv

# ================================================================
#  USER SETTINGS
# ================================================================
FREQ_FILE  = 'results_4x4x4/phonon/Nb_dense_lowfreq.freq'
OUTPUT_TXT = 'DOS_Full_q_real.txt'

CM1_TO_THZ = 0.02998          # cm^-1 → THz conversion
DELTA_W    = 0.05             # THz — Gaussian broadening sigma
MIN_W      = 0.0              # THz — DOS frequency range start
MAX_W      = 7.0              # THz — DOS frequency range end
N_W        = 2000             # number of frequency points

# Branch labels and colors (extend if nbnd > 3)
BRANCH_LABELS = ['TA$_1$', 'TA$_2$', 'LA']
BRANCH_COLORS = ['#2196F3', '#FF5722', '#4CAF50']

# ================================================================
#  PARSE .freq FILE
# ================================================================
def parse_freq_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find &plot header and extract nbnd, nks
    # Handles: single-line, multi-line, and old format
    nbnd, nks, data_start = None, None, 0
    for idx, line in enumerate(lines):
        m_nbnd = re.search(r'nbnd\s*=\s*(\d+)', line, re.IGNORECASE)
        m_nks  = re.search(r'nks\s*=\s*(\d+)',  line, re.IGNORECASE)
        if m_nbnd: nbnd = int(m_nbnd.group(1))
        if m_nks:  nks  = int(m_nks.group(1))
        # Old format: first line is just "nbnd nks"
        if idx == 0 and re.match(r'^\s*\d+\s+\d+\s*$', line):
            nbnd, nks = map(int, line.split())
            data_start = 1
            break
        if '/' in line:
            data_start = idx + 1
            break

    # If nks still None — count q-points from data lines directly
    if nks is None:
        data_lines_tmp = [l.strip() for l in lines[data_start:] if l.strip()]
        # Each q-point = 1 coord line + 1 freq line = 2 lines
        nks = len(data_lines_tmp) // 2
        print(f"  nks not in header — counted {nks} q-points from data")

    print(f"  nbnd={nbnd}, nks={nks}, data starts at line {data_start+1}")

    # Parse q-points and frequencies
    data_lines = [l.strip() for l in lines[data_start:] if l.strip()]

    qpoints = []
    freqs   = []   # shape: (nks, nbnd)

    for k in range(nks):
        coord_line = data_lines[2 * k]
        freq_line  = data_lines[2 * k + 1]
        qpoints.append(tuple(map(float, coord_line.split())))
        freqs.append(list(map(float, freq_line.split())))

    return np.array(qpoints), np.array(freqs), nbnd, nks


# ================================================================
#  GAUSSIAN DOS  (vectorized — fast)
# ================================================================
def gaussian_dos(phonon_freqs_thz, W, sigma):
    """
    phonon_freqs_thz : 1D array of frequencies for one branch (THz)
    W                : frequency axis array (THz)
    sigma            : Gaussian broadening (THz)
    Returns DOS array same length as W
    """
    # Shape: (len(W), len(phonon_freqs))
    diff    = W[:, np.newaxis] - phonon_freqs_thz[np.newaxis, :]
    gauss   = np.exp(-0.5 * (diff / sigma)**2)
    gauss  /= (sigma * np.sqrt(2 * np.pi))
    return gauss.sum(axis=1) / len(phonon_freqs_thz)


# ================================================================
#  MAIN
# ================================================================
print("Loading phonon frequencies...")
qpoints, freqs_cm1, nbnd, nks = parse_freq_file(FREQ_FILE)
print(f"  Loaded {nks} q-points, {nbnd} branches")

# Convert cm^-1 → THz
freqs_thz = freqs_cm1 * CM1_TO_THZ
print(f"  Max frequency: {freqs_thz.max():.3f} THz")
print(f"  Min frequency: {freqs_thz.min():.3f} THz")

# Frequency axis
W = np.linspace(MIN_W, MAX_W, N_W)

# Compute branch DOS
print("Computing branch DOS...")
branch_dos = np.zeros((nbnd, len(W)))
for ib in range(nbnd):
    # Skip negative frequencies (numerical artifacts)
    branch_freqs = freqs_thz[:, ib]
    branch_freqs = branch_freqs[branch_freqs > 0]
    branch_dos[ib] = gaussian_dos(branch_freqs, W, DELTA_W)

total_dos = branch_dos.sum(axis=0)
print("  Done ✅")

# ================================================================
#  NORMALIZED FRACTIONS
# ================================================================
total_safe = np.where(total_dos > 1e-10, total_dos, np.nan)
frac = np.zeros((nbnd, len(W)))
for ib in range(nbnd):
    frac[ib] = branch_dos[ib] / total_safe

# ================================================================
#  PLOT — 3 panels
#  Top    : Full range branch DOS + total
#  Middle : Zoom 0-1 THz
#  Bottom : Normalized fractions 0-1 THz
# ================================================================
fig, axes = plt.subplots(3, 1, figsize=(9, 11))
fig.suptitle("Branch-Resolved Phonon DOS of Nb (BCC)\n"
             "PBE | GBRV-USPP | 4×4×4 q-mesh",
             fontsize=12, fontweight='bold')

# Extend labels/colors if nbnd > 3
labels = BRANCH_LABELS[:nbnd]
colors = BRANCH_COLORS[:nbnd]

# ---- Panel 1: Full range ----
ax = axes[0]
ax.plot(W, total_dos, color='black', linewidth=2.0,
        label='Total', zorder=5)
for ib in range(nbnd):
    ax.plot(W, branch_dos[ib], color=colors[ib],
            linewidth=1.5, linestyle='--',
            alpha=0.9, label=labels[ib])
ax.axvline(x=1.0, color='red', linewidth=1.0,
           linestyle=':', alpha=0.6)
ax.set_xlim(MIN_W, MAX_W)
ax.set_ylim(bottom=0)
ax.set_xlabel("Frequency (THz)", fontsize=11)
ax.set_ylabel("DOS (arb. units)", fontsize=11)
ax.set_title("Full Range", fontsize=11)
ax.legend(fontsize=10, framealpha=0.8)
ax.grid(True, alpha=0.2, linestyle=':')

# ---- Panel 2: Zoom 0-1 THz ----
ax = axes[1]
mask = W <= 1.0
ax.plot(W[mask], total_dos[mask], color='black',
        linewidth=2.0, label='Total')
for ib in range(nbnd):
    ax.plot(W[mask], branch_dos[ib][mask],
            color=colors[ib], linewidth=1.5,
            linestyle='--', alpha=0.9, label=labels[ib])
ax.set_xlim(0, 1.0)
ax.set_ylim(bottom=0)
ax.set_xlabel("Frequency (THz)", fontsize=11)
ax.set_ylabel("DOS (arb. units)", fontsize=11)
ax.set_title("Low Frequency Zoom: 0 – 1 THz", fontsize=11)
ax.legend(fontsize=10, framealpha=0.8)
ax.grid(True, alpha=0.2, linestyle=':')

# ---- Panel 3: Normalized fractions 0-1 THz ----
ax = axes[2]
for ib in range(nbnd):
    ax.plot(W[mask], frac[ib][mask],
            color=colors[ib], linewidth=1.8,
            label=labels[ib])
ax.axhline(y=1/nbnd, color='gray', linewidth=0.8,
           linestyle=':', alpha=0.6,
           label=f'Equal (1/{nbnd})')
ax.axhline(y=0, color='black', linewidth=0.5, alpha=0.3)
ax.axhline(y=1, color='black', linewidth=0.5, alpha=0.3)
ax.set_xlim(0, 1.0)
ax.set_ylim(-0.05, 1.05)
ax.set_xlabel("Frequency (THz)", fontsize=11)
ax.set_ylabel("Fractional contribution", fontsize=11)
ax.set_title(r"Normalized: $DOS_{branch}\,/\,DOS_{total}$  (0–1 THz)",
             fontsize=11)
ax.legend(fontsize=10, framealpha=0.8)
ax.grid(True, alpha=0.2, linestyle=':')

plt.tight_layout()
out_png = 'results_4x4x4/phonon/Nb_branch_dos_clean.png'
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\n  Plot saved: {out_png}")
plt.show()

# ================================================================
#  SAVE TO TEXT FILE
# ================================================================
with open(OUTPUT_TXT, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    # Header
    writer.writerow(['# Freq(THz)'] +
                    [f'Branch_{i+1}' for i in range(nbnd)] +
                    ['Total'])
    # Data
    for i in range(len(W)):
        row = [f'{W[i]:.4f}'] + \
              [f'{branch_dos[ib][i]:.6f}' for ib in range(nbnd)] + \
              [f'{total_dos[i]:.6f}']
        writer.writerow(row)

print(f"  Data saved: {OUTPUT_TXT}")

# ================================================================
#  FRACTIONAL CONTRIBUTIONS TABLE (0-1 THz)
# ================================================================
print("\n  Branch fractional contributions (0–1 THz):")
header = f"  {'Freq(THz)':<12}" + \
         "".join([f'{lb:>10}' for lb in ['TA1','TA2','LA']]) + \
         f"{'Sum':>8}"
print(header)
print("  " + "-" * 48)
check_freqs = np.arange(0.05, 1.05, 0.05)
for f_check in check_freqs:
    idx = np.argmin(np.abs(W - f_check))
    fracs = [frac[ib][idx] for ib in range(nbnd)]
    s = sum(x for x in fracs if not np.isnan(x))
    row = f"  {f_check:<12.2f}" + \
          "".join([f'{x:>10.3f}' if not np.isnan(x) else '       ---'
                   for x in fracs]) + \
          f"{s:>8.3f}"
    print(row)

# ================================================================
#  LINEAR FIT OF FRACTIONAL DOS IN LOW-FREQUENCY REGIME
#  Fits frac[ib](W) = slope * W + intercept  in flat region
#  Evaluates at W = 1.0 THz to get exotic fractional DOS
# ================================================================

# Settings for linear fit
FIT_MIN   = 0.1    # THz — start of flat region for fitting
FIT_MAX   = 1.2    # THz — end of flat region for fitting
EVAL_AT   = 1.0    # THz — evaluate fitted line here

print("\n" + "="*58)
print("  LINEAR FIT OF FRACTIONAL DOS (flat low-freq regime)")
print(f"  Fit range  : {FIT_MIN} – {FIT_MAX} THz")
print(f"  Evaluate at: {EVAL_AT} THz")
print("="*58)

# Mask for fitting region
fit_mask = (W >= FIT_MIN) & (W <= FIT_MAX) & ~np.isnan(total_safe)

W_fit    = W[fit_mask]
fit_results = {}

fig_fit, axes_fit = plt.subplots(1, 2, figsize=(12, 5))
fig_fit.suptitle(
    "Linear Fit of Branch Fractional DOS — Low Frequency Regime\n"
    f"Nb (BCC) | Fit: {FIT_MIN}–{FIT_MAX} THz | Evaluated at {EVAL_AT} THz",
    fontsize=12, fontweight='bold')

# ---- Left panel: fractional DOS + fits ----
ax_l = axes_fit[0]
for ib in range(nbnd):
    frac_fit = frac[ib][fit_mask]

    # Remove NaN points
    valid = ~np.isnan(frac_fit)
    W_v   = W_fit[valid]
    F_v   = frac_fit[valid]

    if len(W_v) < 2:
        print(f"  {BRANCH_LABELS[ib]}: not enough points for fit")
        continue

    # Linear fit: f(W) = slope * W + intercept
    coeffs    = np.polyfit(W_v, F_v, 1)
    slope     = coeffs[0]
    intercept = coeffs[1]
    fit_func  = np.poly1d(coeffs)

    # Evaluate at EVAL_AT
    val_at_eval = fit_func(EVAL_AT)

    # R² goodness of fit
    residuals = F_v - fit_func(W_v)
    ss_res    = np.sum(residuals**2)
    ss_tot    = np.sum((F_v - np.mean(F_v))**2)
    r2        = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    fit_results[ib] = {
        'slope'    : slope,
        'intercept': intercept,
        'val_1THz' : val_at_eval,
        'r2'       : r2
    }

    # Plot data points in fit range
    ax_l.plot(W[fit_mask], frac[ib][fit_mask],
              'o', color=BRANCH_COLORS[ib],
              markersize=2, alpha=0.5)

    # Plot fit line extended to EVAL_AT
    W_ext  = np.linspace(FIT_MIN, EVAL_AT + 0.05, 200)
    ax_l.plot(W_ext, fit_func(W_ext), '-',
              color=BRANCH_COLORS[ib], linewidth=2.0,
              label=f'{BRANCH_LABELS[ib]}  →  {val_at_eval:.4f} at {EVAL_AT} THz')

    # Mark evaluated point
    ax_l.plot(EVAL_AT, val_at_eval, '*',
              color=BRANCH_COLORS[ib], markersize=14,
              zorder=10)

# Reference lines
ax_l.axvline(x=EVAL_AT, color='gray', linewidth=1.0,
             linestyle='--', alpha=0.7, label=f'{EVAL_AT} THz')
ax_l.axhline(y=1/nbnd, color='gray', linewidth=0.8,
             linestyle=':', alpha=0.5, label=f'Equal (1/{nbnd})')
ax_l.set_xlim(FIT_MIN - 0.05, EVAL_AT + 0.1)
ax_l.set_ylim(-0.05, 1.05)
ax_l.set_xlabel("Frequency (THz)", fontsize=11)
ax_l.set_ylabel("Fractional DOS", fontsize=11)
ax_l.set_title(f"Linear Fit  ({FIT_MIN}–{FIT_MAX} THz)", fontsize=11)
ax_l.legend(fontsize=9, framealpha=0.8)
ax_l.grid(True, alpha=0.2, linestyle=':')

# ---- Right panel: bar chart of values at EVAL_AT ----
ax_r = axes_fit[1]
branch_vals  = [fit_results[ib]['val_1THz'] for ib in fit_results]
branch_names = [BRANCH_LABELS[ib].replace('$','').replace('_','')
                for ib in fit_results]
bars = ax_r.bar(branch_names, branch_vals,
                color=[BRANCH_COLORS[ib] for ib in fit_results],
                edgecolor='black', linewidth=0.8, alpha=0.85)

# Annotate bars
for bar, val in zip(bars, branch_vals):
    ax_r.text(bar.get_x() + bar.get_width()/2,
              bar.get_height() + 0.01,
              f'{val:.4f}', ha='center', va='bottom',
              fontsize=11, fontweight='bold')

ax_r.axhline(y=1/nbnd, color='gray', linewidth=1.0,
             linestyle='--', alpha=0.7,
             label=f'Equal contribution (1/{nbnd} = {1/nbnd:.4f})')
ax_r.set_ylim(0, 1.05)
ax_r.set_xlabel("Branch", fontsize=11)
ax_r.set_ylabel(f"Fractional DOS at {EVAL_AT} THz", fontsize=11)
ax_r.set_title(f"Branch Contributions at {EVAL_AT} THz", fontsize=11)
ax_r.legend(fontsize=9, framealpha=0.8)
ax_r.grid(True, axis='y', alpha=0.2, linestyle=':')

plt.tight_layout()
out_fit = 'results_4x4x4/phonon/Nb_branch_dos_fit.png'
plt.savefig(out_fit, dpi=150, bbox_inches='tight')
print(f"  Fit plot saved: {out_fit}")
plt.show()

# ---- Print fit results table ----
print(f"\n  {'Branch':<10} {'Slope':>12} {'Intercept':>12} "
      f"{'Val@{:.1f}THz'.format(EVAL_AT):>14} {'R²':>8}")
print("  " + "-"*58)
for ib in fit_results:
    r = fit_results[ib]
    label = BRANCH_LABELS[ib].replace('$','').replace('_{','').replace('}','')
    print(f"  {label:<10} {r['slope']:>12.6f} {r['intercept']:>12.6f} "
          f"{r['val_1THz']:>14.6f} {r['r2']:>8.6f}")

# Sum check
total_at_eval = sum(fit_results[ib]['val_1THz'] for ib in fit_results)
print(f"  {'Sum':<10} {'---':>12} {'---':>12} "
      f"{total_at_eval:>14.6f} {'---':>8}")
print(f"\n  NOTE: Sum should be ~1.0 — actual: {total_at_eval:.6f}")
