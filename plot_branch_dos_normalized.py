#!/usr/bin/env python3
# ================================================================
#  plot_branch_dos_normalized.py
#  Branch-resolved phonon DOS for Nb (BCC)
#  Top panel    : absolute DOS per branch + total (full range)
#  Bottom panel : normalized branch fraction 0-1 THz
#                 TA1/(total), TA2/(total), LA/(total)
# ================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os, re

RESULTS_DIR = "results_4x4x4/phonon"
FREQ_FILE   = f"{RESULTS_DIR}/Nb.freq"
DOS_FILE    = f"{RESULTS_DIR}/Nb.dos"

CM1_TO_THZ  = 1.0 / 33.3564

BRANCH_LABELS = ['TA$_1$', 'TA$_2$', 'LA']
BRANCH_COLORS = ['#2196F3', '#FF5722', '#4CAF50']   # blue, orange, green
TOTAL_COLOR   = 'black'

# ================================================================
#  Parse .freq file
# ================================================================
def parse_freq_file(filename):
    if not os.path.exists(filename):
        print(f"ERROR: {filename} not found.")
        return None, None, None
    with open(filename) as f:
        lines = f.readlines()
    nbnd, nks, data_start = None, None, 0
    for idx, line in enumerate(lines):
        m_nbnd = re.search(r'nbnd\s*=\s*(\d+)', line, re.IGNORECASE)
        m_nks  = re.search(r'nks\s*=\s*(\d+)',  line, re.IGNORECASE)
        if m_nbnd: nbnd = int(m_nbnd.group(1))
        if m_nks:  nks  = int(m_nks.group(1))
        if idx == 0 and re.match(r'^\s*\d+\s+\d+\s*$', line):
            parts = line.split()
            nbnd, nks = int(parts[0]), int(parts[1])
            data_start = 1
            break
        if '/' in line:
            data_start = idx + 1
            break
    if nbnd is None or nks is None:
        print("ERROR: Could not parse freq header")
        return None, None, None
    q_dist, freqs = [], []
    q_prev, q_cum = None, 0.0
    i = data_start
    for iq in range(nks):
        while i < len(lines) and lines[i].strip() == '':
            i += 1
        if i >= len(lines): break
        try:
            qp = lines[i].split()
            qx, qy, qz = float(qp[0]), float(qp[1]), float(qp[2])
            i += 1
        except (ValueError, IndexError):
            break
        if q_prev is not None:
            dq = np.sqrt((qx-q_prev[0])**2 +
                         (qy-q_prev[1])**2 +
                         (qz-q_prev[2])**2)
            q_cum += dq
        q_dist.append(q_cum)
        q_prev = [qx, qy, qz]
        freq_vals = []
        while len(freq_vals) < nbnd and i < len(lines):
            ln = lines[i].strip()
            i += 1
            if ln == '': continue
            try:
                freq_vals.extend([float(x) for x in ln.split()])
            except ValueError:
                break
        freqs.append(freq_vals[:nbnd])
    return np.array(q_dist), np.array(freqs), nbnd


# ================================================================
#  Parse total DOS file
# ================================================================
def parse_dos_file(filename):
    if not os.path.exists(filename):
        return None, None
    fd, d = [], []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            p = line.split()
            if len(p) >= 2:
                try:
                    fd.append(float(p[0]))
                    d.append(float(p[1]))
                except ValueError:
                    continue
    return np.array(fd), np.array(d)


# ================================================================
#  Compute branch DOS with Gaussian smearing
#  Uses very fine grid + narrow sigma for low-freq resolution
# ================================================================
def compute_branch_dos(freqs_thz, nbnd, freq_min=0.0, freq_max=7.0,
                       deltaE=0.01, sigma=0.05):
    """
    deltaE = 0.01 THz  → very fine grid (700 points over 7 THz)
    sigma  = 0.05 THz  → narrow Gaussian for low-freq detail
    """
    freq_axis  = np.arange(freq_min, freq_max + deltaE, deltaE)
    branch_dos = np.zeros((nbnd, len(freq_axis)))

    for ib in range(nbnd):
        for freq in freqs_thz[:, ib]:
            if freq < 0:
                continue
            gaussian = np.exp(-0.5 * ((freq_axis - freq) / sigma)**2)
            gaussian /= (sigma * np.sqrt(2 * np.pi))
            branch_dos[ib] += gaussian
        branch_dos[ib] /= len(freqs_thz)

    total_dos = branch_dos.sum(axis=0)
    return freq_axis, branch_dos, total_dos


# ================================================================
#  Load
# ================================================================
print("Loading phonon data...")
q_dist, freqs, nbnd = parse_freq_file(FREQ_FILE)
freq_dos_mat, dos_mat = parse_dos_file(DOS_FILE)

if q_dist is None or len(q_dist) == 0:
    print("ERROR: Could not load freq file")
    exit(1)

freqs = freqs * CM1_TO_THZ
if freq_dos_mat is not None:
    freq_dos_mat = freq_dos_mat * CM1_TO_THZ

print(f"  Loaded {len(q_dist)} q-points, {nbnd} branches")
print(f"  Max frequency : {freqs.max():.3f} THz")

# Dual-resolution grid:
# 0.0 - 1.0 THz : ultra-fine 0.002 THz steps, sigma=0.02
# 1.0 - max THz : normal     0.02  THz steps, sigma=0.05
freq_lo, branch_lo, total_lo = compute_branch_dos(
    freqs, nbnd,
    freq_min = 0.0,
    freq_max = 1.0,
    deltaE   = 0.002,
    sigma    = 0.02
)
freq_hi, branch_hi, total_hi = compute_branch_dos(
    freqs, nbnd,
    freq_min = 1.001,
    freq_max = freqs.max() + 0.3,
    deltaE   = 0.02,
    sigma    = 0.05
)
# Stitch grids together
freq_axis  = np.concatenate([freq_lo, freq_hi])
branch_dos = np.concatenate([branch_lo, branch_hi], axis=1)
total_dos  = np.concatenate([total_lo, total_hi])

print("  Branch DOS computed ✅")

# ================================================================
#  Normalized fractions (for 0-1 THz panel)
# ================================================================
# Avoid division by zero
total_safe = np.where(total_dos > 1e-10, total_dos, np.nan)
frac = np.zeros_like(branch_dos)
for ib in range(nbnd):
    frac[ib] = branch_dos[ib] / total_safe

# ================================================================
#  Plot: 3 panels
#  Top    : Full range branch DOS + total
#  Middle : Zoom 0-1 THz branch DOS + total
#  Bottom : Normalized fractions 0-1 THz
# ================================================================
fig = plt.figure(figsize=(9, 12))
gs  = gridspec.GridSpec(3, 1, hspace=0.35,
                        height_ratios=[2, 1.5, 1.5])
ax_full = fig.add_subplot(gs[0])
ax_zoom = fig.add_subplot(gs[1])
ax_frac = fig.add_subplot(gs[2])

fig.suptitle("Branch-Resolved Phonon DOS of Nb (BCC)\n"
             "PBE | GBRV-USPP | 12×12×12 k-mesh | 4×4×4 q-mesh",
             fontsize=12, fontweight='bold')

# ----------------------------------------------------------------
# Panel 1 — Full range DOS
# ----------------------------------------------------------------
ax_full.plot(freq_axis, total_dos,
             color=TOTAL_COLOR, linewidth=2.0,
             label='Total', zorder=5)
for ib in range(nbnd):
    ax_full.plot(freq_axis, branch_dos[ib],
                 color=BRANCH_COLORS[ib], linewidth=1.5,
                 linestyle='--', alpha=0.9,
                 label=BRANCH_LABELS[ib])

ax_full.set_xlim(0, freqs.max() + 0.3)
ax_full.set_ylim(bottom=0)
ax_full.set_xlabel("Frequency (THz)", fontsize=11)
ax_full.set_ylabel("DOS (arb. units)", fontsize=11)
ax_full.set_title("Full Range", fontsize=11)
ax_full.legend(fontsize=10, loc='upper left', framealpha=0.8)
ax_full.axvline(x=1.0, color='red', linewidth=1.0,
                linestyle=':', alpha=0.7, label='1 THz')
ax_full.grid(True, alpha=0.2, linestyle=':')

# ----------------------------------------------------------------
# Panel 2 — Zoom 0 to 1 THz
# ----------------------------------------------------------------
mask = freq_axis <= 1.0
ax_zoom.plot(freq_axis[mask], total_dos[mask],
             color=TOTAL_COLOR, linewidth=2.0, label='Total')
for ib in range(nbnd):
    ax_zoom.plot(freq_axis[mask], branch_dos[ib][mask],
                 color=BRANCH_COLORS[ib], linewidth=1.5,
                 linestyle='--', alpha=0.9,
                 label=BRANCH_LABELS[ib])

ax_zoom.set_xlim(0, 1.0)
ax_zoom.set_ylim(bottom=0)
ax_zoom.set_xlabel("Frequency (THz)", fontsize=11)
ax_zoom.set_ylabel("DOS (arb. units)", fontsize=11)
ax_zoom.set_title("Low Frequency Zoom: 0 – 1 THz", fontsize=11)
ax_zoom.legend(fontsize=10, loc='upper left', framealpha=0.8)
ax_zoom.grid(True, alpha=0.2, linestyle=':')

# ----------------------------------------------------------------
# Panel 3 — Normalized fractions 0 to 1 THz
# ----------------------------------------------------------------
for ib in range(nbnd):
    ax_frac.plot(freq_axis[mask], frac[ib][mask],
                 color=BRANCH_COLORS[ib], linewidth=1.8,
                 label=BRANCH_LABELS[ib])

ax_frac.axhline(y=1/3, color='gray', linewidth=0.8,
                linestyle=':', alpha=0.6, label='1/3')
ax_frac.set_xlim(0, 1.0)
ax_frac.set_ylim(-0.05, 1.05)
ax_frac.set_xlabel("Frequency (THz)", fontsize=11)
ax_frac.set_ylabel("Fractional contribution", fontsize=11)
ax_frac.set_title(
    r"Normalized: $\frac{DOS_{branch}}{DOS_{total}}$  (0 – 1 THz)",
    fontsize=11)
ax_frac.legend(fontsize=10, loc='upper right', framealpha=0.8)
ax_frac.grid(True, alpha=0.2, linestyle=':')

# Add 0% and 100% reference lines
ax_frac.axhline(y=0, color='black', linewidth=0.5, alpha=0.3)
ax_frac.axhline(y=1, color='black', linewidth=0.5, alpha=0.3)

plt.tight_layout()
out_png = f"{RESULTS_DIR}/Nb_branch_dos_normalized.png"
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\n  Plot saved: {out_png}")
plt.show()

# ================================================================
#  Print fractional contributions at key frequencies
# ================================================================
print("\n  Branch fractional contributions in 0-1 THz range:")
print(f"  {'Freq(THz)':<12} {'TA1':>8} {'TA2':>8} {'LA':>8} {'Sum':>8}")
print("  " + "-"*44)
for f_check in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
                    0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]:
    idx = np.argmin(np.abs(freq_axis - f_check))
    fracs = [frac[ib][idx] for ib in range(nbnd)]
    s = sum(f for f in fracs if not np.isnan(f))
    fstr = [f'{f:>8.3f}' if not np.isnan(f) else '     ---'
            for f in fracs]
    print(f"  {f_check:<12.1f} {''.join(fstr)} {s:>8.3f}")
