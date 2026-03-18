#!/bin/bash
# ============================================================
#  CONVERGENCE TEST: ecutwfc for Nb
#  Pseudopotential: nb_pbe_v1.uspp.F.UPF (GBRV, 13 valence e-)
#  Run: bash 01_convergence_ecutwfc.sh
# ============================================================
# === M2 Pro Optimized Settings ===
NCORE=6
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
# =================================


# --- USER SETTINGS ------------------------------------------
PW_EXEC="pw.x"               # or full path: /usr/local/bin/pw.x
PSEUDO_DIR="./pseudo"
TMP_DIR="./tmp"
PSEUDO_FILE="nb_pbe_v1.uspp.F.UPF"
NCORE=4                       # adjust to your Mac's cores
# ------------------------------------------------------------

mkdir -p $TMP_DIR results/ecutwfc

# Cutoff values to test (Ry)
CUTOFFS="40 50 60 70 80"

echo "============================================"
echo " Starting ecutwfc convergence test for Nb"
echo "============================================"
echo "ecut(Ry)   Etot(Ry)" > results/ecutwfc/etot_vs_ecut.dat

for ECUT in $CUTOFFS; do

  ERHO=$(echo "$ECUT * 8" | bc)   # 8x for semi-core USPP
  PREFIX="Nb_ecut${ECUT}"

  echo "  Running ecutwfc = ${ECUT} Ry  (ecutrho = ${ERHO} Ry) ..."

  cat > tmp_scf.in << EOF
&CONTROL
  calculation   = 'scf'
  prefix        = '${PREFIX}'
  outdir        = '${TMP_DIR}'
  pseudo_dir    = '${PSEUDO_DIR}'
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
  degauss       = 0.02
/
&ELECTRONS
  conv_thr      = 1.0d-10
  mixing_beta   = 0.7
/
ATOMIC_SPECIES
  Nb  92.906  ${PSEUDO_FILE}
ATOMIC_POSITIONS {alat}
  Nb  0.0  0.0  0.0
K_POINTS {automatic}
  12 12 12  0 0 0
EOF

  mpirun -np $NCORE $PW_EXEC -npool $NCORE -in tmp_scf.in \
         > results/ecutwfc/${PREFIX}.out 2>&1

  # Extract total energy
  ETOT=$(grep "!    total energy" results/ecutwfc/${PREFIX}.out \
         | tail -1 | awk '{print $5}')

  if [ -z "$ETOT" ]; then
    echo "  WARNING: No energy found for ecut=${ECUT}. Check output."
    ETOT="FAILED"
  fi

  echo "  ecutwfc = ${ECUT} Ry  ->  Etot = ${ETOT} Ry"
  echo "${ECUT}   ${ETOT}" >> results/ecutwfc/etot_vs_ecut.dat

done

rm -f tmp_scf.in

echo ""
echo "============================================"
echo " RESULTS saved to: results/ecutwfc/etot_vs_ecut.dat"
echo "============================================"
cat results/ecutwfc/etot_vs_ecut.dat
echo ""
echo " Choose the ecutwfc where |ΔEtot| < 1 meV/atom (~7.4e-5 Ry)"
