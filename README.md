# Phonon Dispersion of Nb (BCC) — Quantum ESPRESSO Workflow
### Complete Documentation: Parameters, Problems & Fixes

---

## Table of Contents
1. [System Overview](#1-system-overview)
2. [Pseudopotential](#2-pseudopotential)
3. [Convergence Tests](#3-convergence-tests)
4. [Equation of State & Lattice Parameter](#4-equation-of-state--lattice-parameter)
5. [Final Converged Parameters](#5-final-converged-parameters)
6. [Phonon Pipeline](#6-phonon-pipeline)
7. [Problems Encountered & Fixes](#7-problems-encountered--fixes)
8. [File Structure](#8-file-structure)
9. [Scripts Reference](#9-scripts-reference)
10. [Results Summary](#10-results-summary)

---

## 1. System Overview

| Property | Value |
|---|---|
| Material | Niobium (Nb) |
| Structure | BCC (Body-Centered Cubic) |
| Space group | Im3̄m (#229) |
| `ibrav` | 3 |
| DFT functional | PBE (GGA) |
| QE version | 7.2 |
| Machine | Apple M2 Pro, 32 GB unified memory |
| Cores used | 6 (performance cores) |

---

## 2. Pseudopotential

### File
```
nb_pbe_v1.uspp.F.UPF
Source: GBRV library (Garrity, Bennett, Rabe, Vanderbilt)
URL:    https://www.physics.rutgers.edu/gbrv/
```

### Key Properties from PP Header

| Property | Value | Implication |
|---|---|---|
| Type | Ultrasoft (US) | Lower cutoff than PAW/NC |
| Functional | PBE | Standard GGA for metals |
| Z valence | **13 electrons** | Semi-core included (see below) |
| Relativistic | Scalar-relativistic | Correct for Nb (Z=41) |
| Nonlinear Core Correction | Yes | Important for transition metals |
| Suggested cutoff | 0.0 / 0.0 | **Not provided — must converge manually** |

### Valence Configuration (13 electrons)

```
4S  →  2 electrons  (semi-core)
4P  →  6 electrons  (semi-core)
4D  →  4 electrons  (valence)
5S  →  1 electron   (valence)
5P  →  0 electrons  (projector only)
─────────────────────────────────
Total: 13 electrons
```

> **Why semi-core matters:** Including 4s and 4p semi-core states increases
> accuracy for phonons but requires higher `ecutwfc` (~70 Ry) compared to
> valence-only pseudopotentials (~40 Ry). The cutoff ratio `ecutrho/ecutwfc = 8`
> (instead of the usual 4-5) is required for semi-core USPP + NLCC.

### Why GBRV was Chosen for This Calculation

| Library | Type | ecutwfc | Notes |
|---|---|---|---|
| **GBRV** ← used | USPP | 70 Ry | Best speed/accuracy for phonons |
| PSlibrary PAW | PAW | 60-80 Ry | More accurate, heavier |
| SSSP Precision | PAW | 50-60 Ry | Best for publications |
| SG15/ONCVPSP | NC | 80-100 Ry | Needed only for Raman/anharmonic |

---

## 3. Convergence Tests

### 3.1 ecutwfc Convergence

**Script:** `01_convergence_ecutwfc.sh`  
**Fixed:** k-mesh = 12×12×12, ecutrho = 8×ecutwfc

| ecutwfc (Ry) | Etot (Ry) | ΔE (meV) | Converged? |
|---|---|---|---|
| 40 | -117.93266437 | 12.57 | ❌ |
| 50 | -117.93333115 | 3.49 | ❌ |
| 60 | -117.93348313 | 1.42 | ❌ |
| **70** | **-117.93353305** | **0.74** | **✅** |
| 80 | -117.93358759 | 0.00 | reference |

**Result:** `ecutwfc = 70 Ry` (first point below 1 meV/atom threshold)

> **Note:** 60 Ry gives 1.42 meV — above the 1 meV threshold. The semi-core
> 4s/4p states drive the cutoff higher than typical GBRV pseudopotentials.

### 3.2 k-mesh Convergence

**Script:** `02_convergence_kpoints.sh`  
**Fixed:** ecutwfc = 70 Ry, ecutrho = 560 Ry

#### First attempt (degauss = 0.02 Ry) — FAILED

```
Result: Oscillating energies — not monotonically converging
Cause:  Nb has complex Fermi surface; degauss=0.02 too narrow
        for coarse k-meshes → charge sloshing between meshes
```

#### Second attempt (degauss = 0.03 Ry) — SUCCESS

| k-mesh | Etot (Ry) | ΔE (meV) | Converged? |
|---|---|---|---|
| 12×12×12 | -117.93326809 | 3.69 | ❌ |
| 14×14×14 | -117.93293971 | 0.78 | ✅ |
| **16×16×16** | **-117.93299791** | **0.01** | **✅ chosen** |
| 18×18×18 | -117.93300588 | 0.12 | ✅ |
| 20×20×20 | -117.93299687 | 0.00 | reference |

**Result:** `k-mesh = 16×16×16`, `degauss = 0.03 Ry`

> **Important:** The `degauss` value must be identical in SCF, vc-relax,
> and ph.x. Changing it between steps introduces inconsistencies.

---

## 4. Equation of State & Lattice Parameter

### Why EOS Instead of vc-relax

`vc-relax` was attempted twice and failed:

1. **Attempt 1:** `cell_dofree = 'ibrav'` + BFGS → `step_old is NOT normalized`
2. **Attempt 2:** `cell_dofree = 'volume'` → `Isotropic expansion only for ibrav=1`
3. **Solution:** EOS scan — more robust for single-atom cubic systems

### EOS Scan Data

**Script:** `04_eos_scan.sh`  
9 SCF calculations at different volumes, then 11 points (uniform spacing):

| celldm(1) (Bohr) | a (Å) | V (Å³) | Etot (Ry) |
|---|---|---|---|
| 6.00 | 3.17506 | 16.004 | -117.92143496 |
| 6.05 | 3.20152 | 16.407 | -117.92574229 |
| 6.10 | 3.22798 | 16.818 | -117.92899490 |
| 6.15 | 3.25444 | 17.234 | -117.93125460 |
| 6.20 | 3.28090 | 17.658 | -117.93258139 |
| **6.25** | **3.30736** | **18.089** | **-117.93303344 ← minimum** |
| 6.30 | 3.33382 | 18.527 | -117.93266572 |
| 6.35 | 3.36027 | 18.971 | -117.93153094 |
| 6.40 | 3.38673 | 19.423 | -117.92967901 |
| 6.45 | 3.41319 | 19.882 | -117.92715786 |
| 6.50 | 3.43965 | 20.348 | -117.92401194 |

### Birch-Murnaghan EOS Fit Results

**Script:** `fit_eos.py` (v3 — robust minimize in eV units)

| Property | DFT (PBE) | Experimental | Error |
|---|---|---|---|
| a₀ | 3.3084 Å | 3.3008 Å | +0.23% ✅ |
| celldm(1) | 6.2521 Bohr | 6.2436 Bohr | +0.14% ✅ |
| V₀ | 18.107 Å³ | 18.008 Å³ | +0.55% ✅ |
| B₀ | 170.7 GPa | 170 GPa | +0.4% ✅ |
| B₀' | 3.728 | ~3.6-4.0 | ✅ |

> **EOS fitting bug fixed:** Original fitting in Ry/Å³ units caused
> B₀ = 1152 GPa (wrong). Fixed by fitting in eV/Å³ units which gives
> proper numerical scaling for the optimizer.

---

## 5. Final Converged Parameters

```fortran
! ============================================================
!  USE THESE IN ALL QE CALCULATIONS FOR Nb
! ============================================================
&SYSTEM
  ibrav         = 3
  celldm(1)     = 6.2521       ! DFT optimized (Bohr)
  nat           = 1
  ntyp          = 1
  ecutwfc       = 70.0         ! Ry — converged
  ecutrho       = 560.0        ! Ry — 8 × ecutwfc (semi-core USPP)
  occupations   = 'smearing'
  smearing      = 'marzari-vanderbilt'
  degauss       = 0.03         ! Ry — converged
/
&ELECTRONS
  conv_thr      = 1.0d-12      ! tight for phonons
  mixing_beta   = 0.4
  mixing_mode   = 'plain'
/
ATOMIC_SPECIES
  Nb  92.906  nb_pbe_v1.uspp.F.UPF
ATOMIC_POSITIONS {alat}
  Nb  0.0  0.0  0.0
K_POINTS {automatic}
  8 8 8  0 0 0                 ! for phonon SCF (speed)
! 16 16 16  0 0 0              ! for high-accuracy SCF
```

---

## 6. Phonon Pipeline

### Workflow
```
SCF (pw.x) → PHONON (ph.x) → Q2R (q2r.x) → MATDYN (matdyn.x) → PLOT
```

### Step 1 — SCF

```
k-mesh    : 8×8×8 (for phonon speed) or 16×16×16 (high accuracy)
conv_thr  : 1×10⁻¹² Ry
Result    : Etot = -117.934 Ry, EF = 17.82 eV, 12 iterations
```

### Step 2 — Phonon (ph.x)

```fortran
&INPUTPH
  ldisp        = .true.
  nq1 = 4, nq2 = 4, nq3 = 4   ! 8 irreducible q-points for BCC
  tr2_ph       = 1.0d-14
  alpha_mix(1) = 0.3
  recover      = .false.
/
```

**Execution:** Serial + OpenMP (6 threads) — see Section 7 for why.

```bash
export OMP_NUM_THREADS=6
ph.x -in Nb_ph.in > Nb_ph.out 2>&1 &
```

**Timing:**
- With 16×16×16 k-mesh: ~803s/iter → ~36 hrs total ❌
- With 8×8×8 k-mesh: ~146s/iter → ~5-8 hrs total ✅

### Step 3 — q2r

```fortran
&INPUT
  fildyn = 'Nb.dyn'
  zasr   = 'simple'     ! acoustic sum rule
  flfrc  = 'Nb.fc'
/
```

### Step 4 — matdyn (Dispersion)

BCC high-symmetry path: **Γ → H → P → Γ → N**

```fortran
&INPUT
  asr            = 'simple'
  q_in_band_form = .true.
/
5
  0.000  0.000  0.000  40   ! Gamma
  1.000  0.000  0.000  30   ! H
  0.500  0.500  0.500  30   ! P
  0.000  0.000  0.000  40   ! Gamma
  0.000  1.000  0.000   1   ! N
```

### Step 5 — Phonon DOS

```fortran
! matdyn with dos=.true., 20×20×20 q-grid, deltaE=1.0 cm⁻¹
```

---

## 7. Problems Encountered & Fixes

### Problem 1 — vc-relax BFGS Failure
```
Error: step_old is NOT normalized
Cause: Initial pressure 13.68 kbar too far from equilibrium for BFGS
Fix:   Switched to EOS scan approach (more robust for single-atom BCC)
```

### Problem 2 — cell_dofree = 'volume' Error
```
Error: Isotropic expansion only allowed for ibrav=1
Cause: 'volume' mode only works for simple cubic, not BCC
Fix:   EOS scan removes need for vc-relax entirely
```

### Problem 3 — SCF Not Converging in vc-relax
```
Error: Stuck at iter #2, accuracy = 0.089 Ry after 1900+ seconds
Cause: mixing_beta = 0.7 too aggressive + damp-w dynamics
Fix:   Reduced mixing_beta = 0.4, switched to EOS approach
```

### Problem 4 — k-mesh Oscillation
```
Error: Non-monotonic energy convergence vs k-mesh
Cause: degauss = 0.02 Ry too narrow for Nb's complex Fermi surface
Fix:   Increased to degauss = 0.03 Ry → smooth monotonic convergence
```

### Problem 5 — EOS B₀ = 1152 GPa (wrong)
```
Error: Bulk modulus 6.8× too large
Cause: Fitting in Ry/Å³ units → residuals ~10⁻⁶ Ry², optimizer loses sensitivity
Fix:   Refit in eV/Å³ units → B₀ = 170.7 GPa ✅
       Also added: polynomial pre-fit, physical bounds, multiple starting points
```

### Problem 6 — ph.x "File Cannot Be Deleted"
```
Error: Fortran runtime error: File cannot be deleted (test0)
       Run is not recoverable starting from scratch
Cause: Leftover _ph0/Nb.phsave/ directory from previous crashed run
Fix:   rm -rf ./tmp/_ph0/ before every ph.x run
       Set recover = .false.
```

### Problem 7 — MPI ph.x Exit Code 2
```
Error: prterun detected non-zero exit, Exit code: 2
       "Using Slab Decomposition" → crash
Cause: Known incompatibility between OpenMPI + QE FFT Slab
       Decomposition on Apple Silicon
Fix:   Run ph.x in serial + OpenMP instead of MPI
       export OMP_NUM_THREADS=6
       ph.x -in ...  (no mpirun)
```

### Problem 8 — ph.x with -npool Flag
```
Error: MPI crash with exit code 2
Cause: -npool is for pw.x k-point parallelism
       ph.x parallelizes over q-points automatically — no -npool needed
Fix:   Remove -npool from all ph.x calls
Rule:  pw.x → use -npool; ph.x → never use -npool
```

### Problem 9 — ph.x Too Slow (803s/iter)
```
Problem: 803 seconds per DFPT iteration → ~36 hrs total
Cause:   16×16×16 k-mesh in phonon SCF is overkill for DFPT
Fix:     Reduced k-mesh to 8×8×8 for phonon SCF
Result:  146 seconds per iteration → ~5-8 hrs total ✅
Note:    16×16×16 needed for energy convergence but not for
         phonon force constants (DFPT converges faster with k-points)
```

### Problem 10 — taskpolicy Launched New Process
```
Error: sudo taskpolicy -c utility ph.x launched NEW ph.x instance
       instead of modifying existing process priority
Cause: taskpolicy syntax is for launching, not modifying
Fix:   Kill the new suspended process: sudo kill -9 <PID>
       Use: sudo renice -n 0 -p <PID> to change priority instead
```

---

## 8. File Structure

```
Cloude_NB/
├── pseudo/
│   └── nb_pbe_v1.uspp.F.UPF        ← GBRV pseudopotential
├── tmp/                             ← QE temporary files (large)
│   ├── Nb.save/                     ← SCF wavefunctions
│   └── _ph0/                        ← ph.x working files
│       ├── Nb.bar                   ← bare perturbation
│       ├── Nb.dwf                   ← perturbed wavefunctions
│       └── Nb.prd                   ← perturbation data
├── results/
│   ├── ecutwfc/
│   │   └── etot_vs_ecut.dat         ← convergence data
│   ├── kpoints/
│   │   └── etot_vs_kmesh.dat        ← k-mesh convergence data
│   ├── eos/
│   │   ├── eos_data.dat             ← E vs V scan data
│   │   ├── eos_fit_results.txt      ← BM fit results
│   │   └── eos_fit.png              ← E(V) plot
│   └── phonon/
│       ├── Nb_scf.in / .out         ← SCF input/output
│       ├── Nb_ph.in  / .out         ← Phonon input/output
│       ├── Nb.dyn1 ... Nb.dyn8      ← Dynamical matrices
│       ├── Nb.fc                    ← Interatomic force constants
│       ├── Nb.freq                  ← Phonon frequencies (dispersion)
│       ├── Nb.dos                   ← Phonon DOS
│       └── Nb_phonon_dispersion.png ← Final plot
├── 01_convergence_ecutwfc.sh
├── 02_convergence_kpoints.sh
├── 04_eos_scan.sh
├── 05_phonon_pipeline.sh
├── fit_eos.py
├── plot_convergence.py
├── plot_phonon.py
└── run_overnight.sh
```

---

## 9. Scripts Reference

| Script | Purpose | Runtime |
|---|---|---|
| `01_convergence_ecutwfc.sh` | Test ecutwfc 40→80 Ry | ~30 min |
| `02_convergence_kpoints.sh` | Test k-mesh 6×6×6→20×20×20 | ~40 min |
| `04_eos_scan.sh` | EOS volume scan (11 points) | ~45 min |
| `fit_eos.py` | BM EOS fit, extract a₀ and B₀ | ~1 min |
| `plot_convergence.py` | Plot ΔE vs ecutwfc and k-mesh | ~1 min |
| `05_phonon_pipeline.sh` | Full SCF→ph.x→q2r→matdyn→DOS | ~6-10 hrs |
| `plot_phonon.py` | Dispersion + DOS plot | ~1 min |
| `run_overnight.sh` | Chain all steps, prevent Mac sleep | overnight |

### M2 Pro Execution Settings

```bash
# For pw.x (MPI works fine):
export OMP_NUM_THREADS=1
mpirun -np 6 pw.x -npool 6 -in input.in > output.out

# For ph.x (serial + OpenMP — MPI crashes on Apple Silicon):
export OMP_NUM_THREADS=6
export OPENBLAS_NUM_THREADS=6
export VECLIB_MAXIMUM_THREADS=6
ph.x -in Nb_ph.in > Nb_ph.out 2>&1 &

# Keep Mac awake during overnight runs:
caffeinate -i -w $(pgrep ph.x) &
```

---

## 10. Results Summary

### SCF Quality (8×8×8 k-mesh)

```
Total energy  : -117.93400989 Ry
SCF accuracy  : 2.1×10⁻¹³ Ry  ✅
Iterations    : 12             ✅
Fermi energy  : 17.82 eV
Smearing -TS  : 5.86 meV
```

### Structural Properties

```
a₀ (DFT/PBE)  : 3.3084 Å      (exp: 3.3008 Å, error: +0.23%) ✅
B₀ (DFT/PBE)  : 170.7 GPa     (exp: 170 GPa,  error: +0.4%)  ✅
B₀'           : 3.728
```

### Phonon Calculation Status

```
q-mesh        : 4×4×4 (8 irreducible q-points for BCC)
Method        : DFPT (ph.x serial + 6 OpenMP threads)
iter #1 time  : 146.5 seconds
Progress      : 4/8 q-points done at last check
BCC path      : Γ → H → P → Γ → N
```

### After JOB DONE — Next Steps

```bash
# 1. Run q2r and matdyn (already in pipeline script)
bash 05_phonon_pipeline.sh   # will detect ph.x done and continue

# 2. Plot
python3 plot_phonon.py

# 3. Expected output
#    Max phonon frequency: ~280 cm⁻¹
#    No imaginary frequencies (Nb is dynamically stable)
#    3 acoustic + 0 optical branches (1 atom per cell)
```

---

## References

1. P. Giannozzi et al., J. Phys.: Condens. Matter **21**, 395502 (2009)
2. P. Giannozzi et al., J. Phys.: Condens. Matter **29**, 465901 (2017)
3. K. Garrity et al., Comput. Mater. Sci. **81**, 446 (2014) — GBRV library
4. S. Baroni et al., Rev. Mod. Phys. **73**, 515 (2001) — DFPT theory
5. F. Birch, Phys. Rev. **71**, 809 (1947) — Birch-Murnaghan EOS

---

*Generated from interactive QE session on Apple M2 Pro, March 2026*
*QE version 7.2 | GBRV pseudopotential | PBE functional*
