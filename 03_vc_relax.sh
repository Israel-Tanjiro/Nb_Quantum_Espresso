#!/bin/bash
# ============================================================
#  vc-relax: Variable-cell relaxation for Nb
#  Obtains the DFT-optimized lattice parameter
#  Run AFTER scripts 01 and 02
#  Set converged values below before running
#  Run: bash 03_vc_relax.sh
# ============================================================
# === M2 Pro Optimized Settings ===
NCORE=6
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
# =================================
# --- USER SETTINGS ------------------------------------------
PW_EXEC="pw.x"
PSEUDO_DIR="./pseudo"
TMP_DIR="./tmp"
PSEUDO_FILE="nb_pbe_v1.uspp.F.UPF"
# !! SET THESE after running scripts 01 and 02 !!
ECUT=70                       # Ry — your converged ecutwfc
ERHO=560
KMESH=16                       # your converged k-mesh size
# ------------------------------------------------------------

mkdir -p $TMP_DIR results/vcrelax

PREFIX="Nb_vcrelax"

echo "============================================"
echo " Starting vc-relax for Nb"
echo " ecutwfc = ${ECUT} Ry | k-mesh = ${KMESH}x${KMESH}x${KMESH}"
echo "============================================"

cat > results/vcrelax/${PREFIX}.in << EOF
&CONTROL
  calculation   = 'vc-relax'
  prefix        = '${PREFIX}'
  outdir        = '${TMP_DIR}'
  pseudo_dir    = '${PSEUDO_DIR}'
  forc_conv_thr = 1.0d-4
  etot_conv_thr = 1.0d-6
/
&SYSTEM
  ibrav         = 3
  celldm(1)     = 6.235
  nat           = 1
  ntyp          = 1
  ecutwfc       = ${ECUT}
  ecutrho       = ${ERHO}
  occupations   = 'smearing'
  smearing      = 'marzari-vanderbilt'
  degauss       = 0.03
/
&ELECTRONS
  conv_thr      = 1.0d-10
  mixing_beta   = 0.7
/
&IONS
  ion_dynamics  = 'bfgs'
/
&CELL
  cell_dynamics = 'bfgs'
  press_conv_thr = 0.1          ! kbar
  cell_dofree   = 'ibrav'       ! only relax volume, keep BCC symmetry
/
ATOMIC_SPECIES
  Nb  92.906  ${PSEUDO_FILE}
ATOMIC_POSITIONS {alat}
  Nb  0.0  0.0  0.0
K_POINTS {automatic}
  ${KMESH} ${KMESH} ${KMESH}  0 0 0
EOF

echo " Running vc-relax ..."
mpirun -np $NCORE $PW_EXEC -npool $NCORE \
       -in  results/vcrelax/${PREFIX}.in \
       -out results/vcrelax/${PREFIX}.out 2>&1

# ---- Extract results ----------------------------------------
echo ""
echo "============================================"
echo " RESULTS"
echo "============================================"

# Final lattice parameter (celldm(1) in Bohr)
CELLDM=$(grep "celldm(1)" results/vcrelax/${PREFIX}.out \
         | tail -1 | awk '{print $2}')

# Convert to Angstrom (1 Bohr = 0.529177 Ang)
if [ ! -z "$CELLDM" ]; then
  A_ANG=$(echo "scale=6; $CELLDM * 0.529177" | bc)
  echo " Optimized celldm(1) = ${CELLDM} Bohr"
  echo " Optimized lattice a = ${A_ANG} Angstrom"
else
  echo " Could not extract celldm — check results/vcrelax/${PREFIX}.out"
fi

# Final total energy
ETOT=$(grep "!    total energy" results/vcrelax/${PREFIX}.out \
       | tail -1 | awk '{print $5}')
echo " Final total energy   = ${ETOT} Ry"

# Final pressure
PRESS=$(grep "P=" results/vcrelax/${PREFIX}.out \
        | tail -1 | awk '{print $NF}')
echo " Final pressure       = ${PRESS} kbar"

echo ""
echo " Full output: results/vcrelax/${PREFIX}.out"
echo ""
echo " !! Use the optimized celldm(1) in ALL subsequent calculations !!"
echo " !! (SCF, phonon) — do NOT use the experimental value          !!"
