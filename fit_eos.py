#!/usr/bin/env python3
# ============================================================
#  fit_eos.py  (v3 — robust minimize approach)
#  Fits E(V) data to Birch-Murnaghan 3rd order EOS
#  Run: python3 fit_eos.py
#  Requires: numpy, scipy, matplotlib
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, brentq
import os

# Conversion: 1 Ry/Ang^3 = 14710.5 GPa
RY_ANG3_TO_GPA = 14710.5

# ============================================================
#  BM EOS — B0 in GPa for numerical stability
# ============================================================
def birch_murnaghan(V, E0, V0, B0_GPa, Bp):
    B0 = B0_GPa / RY_ANG3_TO_GPA
    eta = (V0 / V) ** (2.0 / 3.0)
    return E0 + (9.0 * V0 * B0 / 16.0) * (
        (eta - 1.0)**3 * Bp +
        (eta - 1.0)**2 * (6.0 - 4.0 * eta)
    )

def residuals_sq(params, V, E):
    return np.sum((E - birch_murnaghan(V, *params))**2)

# ============================================================
#  Load data
# ============================================================
data_file = "results/eos/eos_data.dat"

if not os.path.exists(data_file):
    print(f"ERROR: {data_file} not found. Run 04_eos_scan.sh first.")
    exit(1)

volumes, energies, lattices = [], [], []

with open(data_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) == 4 and parts[3] != "FAILED":
            lattices.append(float(parts[1]))
            volumes.append(float(parts[2]))
            energies.append(float(parts[3]))

V = np.array(volumes)
E = np.array(energies)
a = np.array(lattices)

print(f"Loaded {len(V)} data points")
print(f"Volume range : {V.min():.3f} — {V.max():.3f} Ang^3")
print(f"Energy range : {E.min():.8f} — {E.max():.8f} Ry")
print()

# ============================================================
#  Step 1 — Polynomial pre-fit for robust initial guess
#  B0 = V0 * d2E/dV2
# ============================================================
poly   = np.polyfit(V, E, 4)
d1poly = np.polyder(poly, 1)
d2poly = np.polyder(poly, 2)

V0_poly = brentq(lambda v: np.polyval(d1poly, v), V.min(), V.max())
E0_poly = np.polyval(poly, V0_poly)
d2E     = np.polyval(d2poly, V0_poly)
B0_poly_GPa = V0_poly * d2E * RY_ANG3_TO_GPA

print(f"Polynomial pre-fit (numerical):")
print(f"  V0 = {V0_poly:.4f} Ang^3")
print(f"  E0 = {E0_poly:.8f} Ry")
print(f"  B0 = {B0_poly_GPa:.1f} GPa")
print()

# ============================================================
#  Step 2 — BM fit with multiple starting points
# ============================================================
best_result   = None
best_residual = np.inf

B0_starts = [B0_poly_GPa, 150.0, 170.0, 190.0, 210.0]
Bp_starts = [3.5, 3.8, 4.0, 4.5]

for B0_s in B0_starts:
    for Bp_s in Bp_starts:
        p0 = [E0_poly, V0_poly, B0_s, Bp_s]
        try:
            result = minimize(
                residuals_sq, p0, args=(V, E),
                method='Nelder-Mead',
                options={
                    'maxiter': 100000,
                    'xatol':   1e-12,
                    'fatol':   1e-20,
                    'adaptive': True
                }
            )
            E0_t, V0_t, B0_t, Bp_t = result.x
            physical = (50 < B0_t < 400) and (1.5 < Bp_t < 8.0) \
                       and (V.min() < V0_t < V.max())
            if result.fun < best_residual and physical:
                best_residual = result.fun
                best_result   = result
        except Exception:
            continue

# ============================================================
#  Results
# ============================================================
print("=" * 52)

if best_result is not None:
    E0_fit, V0_fit, B0_fit_GPa, Bp_fit = best_result.x
    a0_ang  = (2.0 * V0_fit) ** (1.0 / 3.0)
    a0_bohr = a0_ang / 0.529177
    E_pred  = birch_murnaghan(V, E0_fit, V0_fit, B0_fit_GPa, Bp_fit)
    rms_err = np.sqrt(np.mean((E - E_pred)**2)) * 1000  # mRy

    print("  BIRCH-MURNAGHAN EOS FIT RESULTS        ")
    print("=" * 52)
    print(f"  E0        = {E0_fit:.8f}  Ry")
    print(f"  V0        = {V0_fit:.5f}    Ang^3")
    print(f"  a0        = {a0_ang:.6f}  Angstrom")
    print(f"  a0        = {a0_bohr:.6f}  Bohr  <- use as celldm(1)")
    print(f"  B0 (BM)   = {B0_fit_GPa:.2f}      GPa")
    print(f"  B0 (poly) = {B0_poly_GPa:.2f}      GPa  (cross-check)")
    print(f"  B0'       = {Bp_fit:.4f}")
    print(f"  RMS error = {rms_err:.5f} mRy")
    print("=" * 52)
    print()
    print(f"  Experimental a0 = 3.3008 Ang")
    print(f"  Your DFT a0     = {a0_ang:.4f} Ang")
    print(f"  Difference      = {abs(a0_ang-3.3008)/3.3008*100:.2f} %")
    print()
    print(f"  Experimental B0 = 170 GPa")
    print(f"  BM  fit  B0     = {B0_fit_GPa:.1f} GPa")
    print(f"  Numerical B0    = {B0_poly_GPa:.1f} GPa")
    print("=" * 52)

    os.makedirs("results/eos", exist_ok=True)
    with open("results/eos/eos_fit_results.txt", "w") as f:
        f.write("BIRCH-MURNAGHAN EOS FIT RESULTS FOR Nb\n")
        f.write("=" * 52 + "\n")
        f.write(f"E0        = {E0_fit:.8f} Ry\n")
        f.write(f"V0        = {V0_fit:.5f} Ang^3\n")
        f.write(f"a0        = {a0_ang:.6f} Angstrom\n")
        f.write(f"celldm(1) = {a0_bohr:.6f} Bohr\n")
        f.write(f"B0 (BM)   = {B0_fit_GPa:.2f} GPa\n")
        f.write(f"B0 (poly) = {B0_poly_GPa:.2f} GPa\n")
        f.write(f"B0'       = {Bp_fit:.4f}\n")
        f.write("=" * 52 + "\n")
        f.write(f"\nUSE IN ALL QE INPUTS:\n")
        f.write(f"  celldm(1) = {a0_bohr:.5f}\n")
    print(f"\n  Saved to: results/eos/eos_fit_results.txt")

else:
    # Fallback to polynomial result
    print("  BM fit did not converge — using polynomial estimate")
    print("=" * 52)
    a0_ang  = (2.0 * V0_poly) ** (1.0 / 3.0)
    a0_bohr = a0_ang / 0.529177
    E0_fit, V0_fit, B0_fit_GPa, Bp_fit = E0_poly, V0_poly, B0_poly_GPa, 3.8
    print(f"  a0        = {a0_ang:.6f}  Angstrom")
    print(f"  a0        = {a0_bohr:.6f}  Bohr  <- use as celldm(1)")
    print(f"  B0 (poly) = {B0_poly_GPa:.2f}      GPa")
    print("=" * 52)

# ============================================================
#  Plot
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle("Equation of State — Nb (BCC, PBE)", fontsize=14, fontweight='bold')

V_fine = np.linspace(V.min()*0.98, V.max()*1.02, 500)
E_bm   = birch_murnaghan(V_fine, E0_fit, V0_fit, B0_fit_GPa, Bp_fit)
a_fine = (2.0 * V_fine) ** (1.0/3.0)

ax = axes[0]
ax.scatter(V, E, color='steelblue', s=90, zorder=5,
           label='DFT data', edgecolors='navy', linewidths=1.5)
ax.plot(V_fine, E_bm, 'r-', linewidth=2, label='BM fit')
ax.axvline(V0_fit, color='green', linestyle='--', linewidth=1.5,
           label=f'V0 = {V0_fit:.3f} A3')
ax.set_xlabel("Volume (A^3/atom)", fontsize=12)
ax.set_ylabel("Total Energy (Ry)", fontsize=12)
ax.set_title("E vs V", fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.scatter(a, E, color='darkorange', s=90, zorder=5,
           label='DFT data', edgecolors='saddlebrown', linewidths=1.5)
ax.plot(a_fine, E_bm, 'r-', linewidth=2, label='BM fit')
ax.axvline(a0_ang, color='green', linestyle='--', linewidth=1.5,
           label=f'a0 = {a0_ang:.4f} A (DFT)')
ax.axvline(3.3008, color='purple', linestyle=':', linewidth=1.5,
           label='a0 = 3.3008 A (Exp.)')
ax.set_xlabel("Lattice parameter a (A)", fontsize=12)
ax.set_ylabel("Total Energy (Ry)", fontsize=12)
ax.set_title("E vs a", fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs("results/eos", exist_ok=True)
plt.savefig("results/eos/eos_fit.png", dpi=150, bbox_inches='tight')
print(f"\n  Plot saved to: results/eos/eos_fit.png")
plt.show()
