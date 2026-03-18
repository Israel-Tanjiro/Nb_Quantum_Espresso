#!/usr/bin/env python3
# ================================================================
#  plot_phonon_branch_dos.py
#  Phonon dispersion + branch-resolved DOS for Nb (BCC)
#  Computes branch DOS directly from Nb.freq file
#  Units: THz
# ================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os, re

RESULTS_DIR = "results_4x4x4/phonon"
FREQ_FILE   = f"{RESULTS_DIR}/Nb.freq"
DOS_FILE    = f"{RESULTS_DIR}/Nb.dos"   # total DOS from matdyn

CM1_TO_THZ  = 1.0 / 33.3564

# Branch labels and colors for Nb BCC (3 acoustic branches)
BRANCH_LABELS = ['TA$_1$', 'TA$_2$', 'LA']
BRANCH_COLORS = ['#2196F3', '#FF5722', '#4CAF50']   # blue, orange, green
TOTAL_COLOR   = 'black'

# ================================================================
#  Parse .freq file (handles &plot single-line header)
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
#  Parse total DOS from matdyn .dos file
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
#  Compute branch-resolved DOS by histogramming each branch
# ================================================================
def compute_branch_dos(freqs_thz, nbnd, freq_min, freq_max,
                       deltaE=0.05, sigma=0.08):
    """
    Gaussian-smeared histogram of each phonon branch.
    freqs_thz : (nq, nbnd) array in THz
    sigma      : Gaussian broadening in THz
    deltaE     : frequency bin width in THz
    """
    freq_axis = np.arange(freq_min, freq_max + deltaE, deltaE)
    branch_dos = np.zeros((nbnd, len(freq_axis)))

    for ib in range(nbnd):
        for freq in freqs_thz[:, ib]:
            if freq < 0:
                continue
            # Gaussian broadening
            gaussian = np.exp(-0.5 * ((freq_axis - freq) / sigma)**2)
            gaussian /= (sigma * np.sqrt(2 * np.pi))
            branch_dos[ib] += gaussian

        # Normalize by number of q-points
        branch_dos[ib] /= len(freqs_thz)

    total_dos = branch_dos.sum(axis=0)
    return freq_axis, branch_dos, total_dos


# ================================================================
#  Load data
# ================================================================
print("Loading phonon data...")
q_dist, freqs, nbnd = parse_freq_file(FREQ_FILE)
freq_dos_mat, dos_mat = parse_dos_file(DOS_FILE)

if q_dist is None or len(q_dist) == 0:
    print("ERROR: Could not load freq file")
    exit(1)

# Convert to THz
freqs = freqs * CM1_TO_THZ
if freq_dos_mat is not None:
    freq_dos_mat = freq_dos_mat * CM1_TO_THZ

print(f"  Loaded {len(q_dist)} q-points, {nbnd} branches")
print(f"  Max frequency : {freqs.max():.3f} THz")

# Floor negatives for display
freqs_plot = np.where(freqs < 0, 0, freqs)

# Compute branch DOS
freq_min  = 0.0
freq_max  = freqs.max() + 0.3
freq_axis, branch_dos, total_dos_calc = compute_branch_dos(
    freqs, nbnd, freq_min, freq_max, deltaE=0.05, sigma=0.08)

print(f"  Branch DOS computed ✅")

# ================================================================
#  High-symmetry points
#  Γ(40) → H(30) → P(30) → Γ(40) → N(1) = 141 points
# ================================================================
seg_counts = [40, 30, 30, 40, 1]
seg_ends   = [0] + list(np.cumsum(seg_counts) - 1)
labels_raw = ['Γ', 'H', 'P', 'Γ', 'N']

hsym_pos, hsym_labels = [], []
for idx, lab in zip(seg_ends, labels_raw):
    if idx < len(q_dist):
        hsym_pos.append(q_dist[idx])
        hsym_labels.append(lab)

# ================================================================
#  Plot: dispersion (left) + branch DOS (right)
# ================================================================
fig = plt.figure(figsize=(12, 6))
gs  = gridspec.GridSpec(1, 2, width_ratios=[3, 1.2], wspace=0.04)
ax_disp = fig.add_subplot(gs[0])
ax_dos  = fig.add_subplot(gs[1])

fig.suptitle(
    "Phonon Dispersion & Branch-Resolved DOS of Nb (BCC)\n"
    "PBE | GBRV-USPP | 12×12×12 k-mesh | 4×4×4 q-mesh",
    fontsize=12, fontweight='bold')

# ---- Dispersion ----
for ib in range(nbnd):
    ax_disp.plot(q_dist, freqs_plot[:, ib], '-',
                 color=BRANCH_COLORS[ib],
                 linewidth=1.8, alpha=0.9,
                 label=BRANCH_LABELS[ib])

for pos in hsym_pos:
    ax_disp.axvline(x=pos, color='gray', linewidth=0.8,
                    linestyle='--', alpha=0.5)
ax_disp.axhline(y=0, color='black', linewidth=0.5, alpha=0.4)

ymax = freqs_plot.max() + 0.5
ax_disp.set_ylim(bottom=-0.2, top=ymax)
ax_disp.set_xlim(q_dist[0], q_dist[-1])
ax_disp.set_xticks(hsym_pos)
ax_disp.set_xticklabels(hsym_labels, fontsize=14)
ax_disp.set_ylabel("Frequency (THz)", fontsize=12)
ax_disp.set_title("Phonon Dispersion", fontsize=11)
ax_disp.grid(True, axis='y', alpha=0.2, linestyle=':')
ax_disp.legend(loc='upper right', fontsize=10, framealpha=0.8)

# Annotate max frequency
imax = np.unravel_index(np.argmax(freqs_plot), freqs_plot.shape)
ax_disp.annotate(f'{freqs_plot.max():.2f} THz',
                 xy=(q_dist[imax[0]], freqs_plot.max()),
                 xytext=(6, 4), textcoords='offset points',
                 fontsize=9, color='navy')

# ---- Branch DOS ----
# Total DOS (from matdyn if available, else computed)
if freq_dos_mat is not None and len(freq_dos_mat) > 0:
    ax_dos.plot(dos_mat, freq_dos_mat,
                color=TOTAL_COLOR, linewidth=1.8,
                label='Total', zorder=5)
else:
    ax_dos.plot(total_dos_calc, freq_axis,
                color=TOTAL_COLOR, linewidth=1.8,
                label='Total', zorder=5)

# Individual branch contributions
for ib in range(nbnd):
    ax_dos.plot(branch_dos[ib], freq_axis,
                color=BRANCH_COLORS[ib],
                linewidth=1.4, linestyle='--',
                alpha=0.9, label=BRANCH_LABELS[ib])

ax_dos.axhline(y=0, color='black', linewidth=0.5, alpha=0.4)
ax_dos.set_ylim(-0.2, ymax)
ax_dos.set_xlim(left=0)
ax_dos.set_xlabel("DOS (arb. units)", fontsize=11)
ax_dos.set_title("Phonon DOS", fontsize=11)
ax_dos.yaxis.set_ticklabels([])
ax_dos.grid(True, axis='y', alpha=0.2, linestyle=':')
ax_dos.legend(loc='upper right', fontsize=9, framealpha=0.8)

plt.tight_layout()
out_png = f"{RESULTS_DIR}/Nb_phonon_branch_dos.png"
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\n  Plot saved: {out_png}")
plt.show()

# ================================================================
#  Branch DOS summary
# ================================================================
print("\n  Branch DOS peak frequencies (THz):")
for ib in range(nbnd):
    peak_idx = np.argmax(branch_dos[ib])
    print(f"  {BRANCH_LABELS[ib]:8s} peak : {freq_axis[peak_idx]:.3f} THz")
print(f"  {'Total':8s} peak : {freq_axis[np.argmax(total_dos_calc)]:.3f} THz")
