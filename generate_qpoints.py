#!/usr/bin/env python3
# ================================================================
#  generate_qpoints.py
#  Generates a q-point mesh for matdyn focused on low frequencies
#
#  Strategy: dense mesh near Gamma (0 to Q_MAX in x,y,z)
#  Output: matdyn input file ready to run directly
#
#  Usage: python3 generate_qpoints.py
#         Then: matdyn.x < Nb_matdyn_dense.in > Nb_matdyn_dense.out
# ================================================================

import numpy as np

# ================================================================
#  USER SETTINGS
# ================================================================
FC_FILE    = 'results_4x4x4/phonon/Nb.fc'
FREQ_FILE  = 'results_4x4x4/phonon/Nb_dense.freq'
VEC_FILE   = 'results_4x4x4/phonon/Nb_dense.modes'
ASR        = 'crystal'

# Q-point mesh settings
Q_MAX      = 0.5    # max q in each direction (0.5 = full BZ)
Q_LOW_MAX  = 0.3    # dense region boundary
N_LOW      = 30     # points from 0 to Q_LOW_MAX (dense)
N_HIGH     = 10     # points from Q_LOW_MAX to Q_MAX (coarse)
WEIGHT     = 1      # weight for all q-points (matdyn convention)

# Output files
QPOINTS_FILE  = 'results_4x4x4/phonon/qpoints_dense.dat'
MATDYN_IN     = 'results_4x4x4/phonon/Nb_matdyn_dense.in'

# ================================================================
#  GENERATE Q-POINT GRID
# ================================================================
def generate_qgrid(q_low_max, n_low, q_max, n_high):
    """
    Generate 1D array of q values:
    - Dense from 0 to q_low_max (n_low points)
    - Coarse from q_low_max to q_max (n_high points)
    """
    q_dense  = np.linspace(0.0,       q_low_max, n_low,  endpoint=False)
    q_coarse = np.linspace(q_low_max, q_max,     n_high, endpoint=True)
    return np.concatenate([q_dense, q_coarse])

# 1D q values
q_vals = generate_qgrid(Q_LOW_MAX, N_LOW, Q_MAX, N_HIGH)
nq_1d  = len(q_vals)

print(f"  1D q-points     : {nq_1d}")
print(f"  Dense region    : {N_LOW} points from 0 to {Q_LOW_MAX}")
print(f"  Coarse region   : {N_HIGH} points from {Q_LOW_MAX} to {Q_MAX}")

# 3D mesh: all combinations of (qx, qy, qz)
qx, qy, qz = np.meshgrid(q_vals, q_vals, q_vals, indexing='ij')
qx = qx.flatten()
qy = qy.flatten()
qz = qz.flatten()

# Remove duplicates due to floating point
qpoints = np.unique(
    np.round(np.column_stack([qx, qy, qz]), 8),
    axis=0
)

nq_total = len(qpoints)
print(f"  3D q-points     : {nq_total}")
print(f"  Low-freq region : q < {Q_LOW_MAX} in all directions")

# ================================================================
#  SAVE Q-POINTS FILE (for reference)
# ================================================================
with open(QPOINTS_FILE, 'w') as f:
    f.write(f"# Q-point mesh: {nq_total} points\n")
    f.write(f"# Dense 0-{Q_LOW_MAX} ({N_LOW} pts), "
            f"Coarse {Q_LOW_MAX}-{Q_MAX} ({N_HIGH} pts)\n")
    for q in qpoints:
        f.write(f"  {q[0]:.10f}  {q[1]:.10f}  {q[2]:.10f}\n")

print(f"  Q-points saved  : {QPOINTS_FILE}")

# ================================================================
#  WRITE MATDYN INPUT FILE
# ================================================================
with open(MATDYN_IN, 'w') as f:
    f.write(f" &INPUT\n")
    f.write(f"  amass(1)         = 92.906\n")
    f.write(f"  flfrc            = '{FC_FILE}'\n")
    f.write(f"  flfrq            = '{FREQ_FILE}'\n")
    f.write(f"  flvec            = '{VEC_FILE}'\n")
    f.write(f"  asr              = '{ASR}'\n")
    f.write(f"  q_in_band_form   = .false.\n")
    f.write(f"  q_in_cryst_coord = .false.\n")
    f.write(f" /\n")
    f.write(f"{nq_total}\n")
    for q in qpoints:
        f.write(f"  {q[0]:.10f}  {q[1]:.10f}  {q[2]:.10f}  {WEIGHT}\n")

print(f"  matdyn input    : {MATDYN_IN}")

# ================================================================
#  SUMMARY
# ================================================================
print(f"\n  Run matdyn with:")
print(f"  matdyn.x < {MATDYN_IN} > results_4x4x4/phonon/Nb_matdyn_dense.out")
print(f"\n  Then compute DOS with:")
print(f"  python3 branch_dos.py  (set FREQ_FILE = '{FREQ_FILE}')")

# ================================================================
#  OPTIONAL: also generate a low-freq ONLY mesh (0 to Q_LOW_MAX)
# ================================================================
mask_low = np.all(qpoints <= Q_LOW_MAX, axis=1)
qpoints_low = qpoints[mask_low]
nq_low = len(qpoints_low)

MATDYN_LOW = MATDYN_IN.replace('.in', '_lowfreq.in')
FREQ_LOW   = FREQ_FILE.replace('.freq', '_lowfreq.freq')
VEC_LOW    = VEC_FILE.replace('.modes', '_lowfreq.modes')

with open(MATDYN_LOW, 'w') as f:
    f.write(f" &INPUT\n")
    f.write(f"  amass(1)         = 92.906\n")
    f.write(f"  flfrc            = '{FC_FILE}'\n")
    f.write(f"  flfrq            = '{FREQ_LOW}'\n")
    f.write(f"  flvec            = '{VEC_LOW}'\n")
    f.write(f"  asr              = '{ASR}'\n")
    f.write(f"  q_in_band_form   = .false.\n")
    f.write(f"  q_in_cryst_coord = .false.\n")
    f.write(f" /\n")
    f.write(f"{nq_low}\n")
    for q in qpoints_low:
        f.write(f"  {q[0]:.10f}  {q[1]:.10f}  {q[2]:.10f}  {WEIGHT}\n")

print(f"\n  Low-freq only mesh ({nq_low} points, q < {Q_LOW_MAX}):")
print(f"  matdyn.x < {MATDYN_LOW} > results_4x4x4/phonon/Nb_matdyn_dense_lowfreq.out")
