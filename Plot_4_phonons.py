#!/usr/bin/env python3
# ================================================================
#  plot_phonon.py  (v5 — fixed for single-line &plot header)
#  Format: "&plot nbnd=   3, nks= 141 /"  on line 1
#          then q-point + freq data from line 2 onwards
# ================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os, re

RESULTS_DIR = "results_4x4x4/phonon"
FREQ_FILE   = f"{RESULTS_DIR}/Nb.freq"
DOS_FILE    = f"{RESULTS_DIR}/Nb.dos"

def parse_freq_file(filename):
    if not os.path.exists(filename):
        print(f"ERROR: {filename} not found.")
        return None, None, None

    with open(filename) as f:
        lines = f.readlines()

    nbnd, nks, data_start = None, None, 0

    for idx, line in enumerate(lines):
        # Try to extract nbnd and nks from this line
        m_nbnd = re.search(r'nbnd\s*=\s*(\d+)', line, re.IGNORECASE)
        m_nks  = re.search(r'nks\s*=\s*(\d+)',  line, re.IGNORECASE)
        if m_nbnd: nbnd = int(m_nbnd.group(1))
        if m_nks:  nks  = int(m_nks.group(1))

        # Old format: "3 141" on first line
        if idx == 0 and re.match(r'^\s*\d+\s+\d+\s*$', line):
            parts = line.split()
            nbnd, nks = int(parts[0]), int(parts[1])
            data_start = 1
            break

        # End of namelist header — check if '/' is on THIS line
        # If so, data starts on the NEXT line
        if '/' in line:
            data_start = idx + 1
            break

    print(f"  nbnd={nbnd}, nks={nks}, data starts at line {data_start+1}")

    if nbnd is None or nks is None:
        print("ERROR: Could not parse header")
        return None, None, None

    # Parse q-points and frequencies
    q_dist, freqs = [], []
    q_prev, q_cum = None, 0.0
    i = data_start

    for iq in range(nks):
        while i < len(lines) and lines[i].strip() == '':
            i += 1
        if i >= len(lines):
            break

        try:
            qp = lines[i].split()
            qx, qy, qz = float(qp[0]), float(qp[1]), float(qp[2])
            i += 1
        except (ValueError, IndexError):
            print(f"  ERROR parsing q at line {i+1}: '{lines[i].strip()}'")
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

    print(f"  Successfully parsed {len(q_dist)} q-points")
    return np.array(q_dist), np.array(freqs), nbnd


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
#  Load
# ================================================================
print("Loading phonon data...")
q_dist, freqs, nbnd = parse_freq_file(FREQ_FILE)
freq_dos, dos = parse_dos_file(DOS_FILE)

if q_dist is None or len(q_dist) == 0:
    print("Failed — run:  cat -n results/phonon/Nb.freq | head -5")
    exit(1)

# Convert cm^-1 → THz  (1 THz = 33.3564 cm^-1)
CM1_TO_THZ = 1.0 / 33.3564
freqs    = freqs    * CM1_TO_THZ
if freq_dos is not None and len(freq_dos) > 0:
    freq_dos = freq_dos * CM1_TO_THZ

print(f"  Max frequency : {freqs.max():.3f} THz")
print(f"  Min frequency : {freqs.min():.3f} THz")
n_imag = np.sum(freqs < -0.15)   # -5 cm^-1 in THz
if n_imag > 0:
    print(f"  WARNING: {n_imag} imaginary frequencies (< -5 cm^-1)")
    print(f"  Setting negative frequencies to 0 for display")
    freqs_plot = np.where(freqs < 0, 0, freqs)
else:
    print(f"  No imaginary frequencies ✅")
    freqs_plot = freqs

# ================================================================
#  High-symmetry points
#  Γ(40) → H(30) → P(30) → Γ(40) → N(1) = 141 points
# ================================================================
seg_counts = [40, 30, 30, 40, 1]
seg_ends   = [0] + list(np.cumsum(seg_counts) - 1)
labels_raw = ['Γ', 'H', 'P', 'Γ', 'N']

hsym_pos, hsym_labels = [], []
seen = []
for idx, lab in zip(seg_ends, labels_raw):
    if idx < len(q_dist):
        pos = q_dist[idx]
        hsym_pos.append(pos)
        hsym_labels.append(lab)

# ================================================================
#  Plot
# ================================================================
has_dos = freq_dos is not None and len(freq_dos) > 0

if has_dos:
    fig = plt.figure(figsize=(10, 6))
    gs  = gridspec.GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05)
    ax_disp = fig.add_subplot(gs[0])
    ax_dos  = fig.add_subplot(gs[1])
else:
    fig, ax_disp = plt.subplots(figsize=(8, 6))

fig.suptitle("Phonon Dispersion of Nb (BCC)\nPBE | GBRV-USPP | 12×12×12 k-mesh | 4×4×4 q-mesh",
             fontsize=12, fontweight='bold')

ax = ax_disp
for ib in range(nbnd):
    ax.plot(q_dist, freqs_plot[:, ib], '-',
            color='steelblue', linewidth=1.8, alpha=0.9)

for pos in hsym_pos:
    ax.axvline(x=pos, color='gray', linewidth=0.8,
               linestyle='--', alpha=0.5)
ax.axhline(y=0, color='black', linewidth=0.6, alpha=0.4,
           linestyle='-')

ax.set_ylim(bottom=-0.2, top=freqs_plot.max()+0.5)
ax.set_xlim(q_dist[0], q_dist[-1])
ax.set_xticks(hsym_pos)
ax.set_xticklabels(hsym_labels, fontsize=14)
ax.set_ylabel("Frequency (THz)", fontsize=12)
ax.set_title("Phonon Dispersion", fontsize=11)
ax.grid(True, axis='y', alpha=0.2, linestyle=':')

imax = np.unravel_index(np.argmax(freqs_plot), freqs_plot.shape)
ax.annotate(f'{freqs_plot.max():.2f} THz',
            xy=(q_dist[imax[0]], freqs_plot.max()),
            xytext=(6, 4), textcoords='offset points',
            fontsize=9, color='navy')

if has_dos:
    ax_dos.plot(dos, freq_dos, 'r-', linewidth=1.5)
    ax_dos.fill_betweenx(freq_dos, 0, dos, alpha=0.25, color='red')
    ax_dos.axhline(y=0, color='black', linewidth=0.6, alpha=0.4)
    ax_dos.set_ylim(ax_disp.get_ylim())
    ax_dos.set_xlabel("DOS", fontsize=11)
    ax_dos.set_title("Phonon DOS", fontsize=11)
    ax_dos.yaxis.set_ticklabels([])
    ax_dos.set_xlim(left=0)
    ax_dos.grid(True, axis='y', alpha=0.2, linestyle=':')

plt.tight_layout()
out_png = f"{RESULTS_DIR}/Nb_phonon_dispersion.png"
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\n  Plot saved: {out_png}")
plt.show()

# Frequency table
print("\n  Frequencies at high-symmetry points (THz):")
print(f"  {'Point':<6} {'Mode 1':>10} {'Mode 2':>10} {'Mode 3':>10}")
print("  " + "-"*38)
for idx, lab in zip(seg_ends, labels_raw):
    if idx < len(q_dist):
        f = freqs[idx]
        print(f"  {lab:<6} {f[0]:>10.3f} {f[1]:>10.3f} {f[2]:>10.3f}")
