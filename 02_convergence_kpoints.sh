#!/bin/bash
# ============================================================
#  CONVERGENCE TEST: k-mesh for Nb
#  Run AFTER 01_convergence_ecutwfc.sh
#  Set ECUT_CONVERGED below before running
#  Run: bash 02_convergence_kpoints.sh
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
NCORE=4

# !! SET THIS after running script 01 !!
ECUT_CONVERGED=70             # Ry — update with your converged value
ERHO= 560
# ------------------------------------------------------------

mkdir -p $TMP_DIR results/kpoints

# k-mesh sizes to test
KMESHES="12 14 16 18 20"

echo "============================================"
echo " Starting k-mesh convergence test for Nb"
echo " Using ecutwfc = ${ECUT_CONVERGED} Ry"
echo "============================================"
echo "k-mesh   Etot(Ry)" > results/kpoints/etot_vs_kmesh.dat

for K in $KMESHES; do

  PREFIX="Nb_k${K}"
  echo "  Running k-mesh = ${K}x${K}x${K} ..."

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
  ecutwfc       = ${ECUT_CONVERGED}
  ecutrho       = ${ERHO}
  occupations   = 'smearing'
  smearing      = 'marzari-vanderbilt'
  degauss       = 0.03
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
  ${K} ${K} ${K}  0 0 0
EOF

  mpirun -np $NCORE $PW_EXEC -npool $NCORE -in tmp_scf.in \
         > results/kpoints/${PREFIX}.out 2>&1

  # Extract total energy
  ETOT=$(grep "!    total energy" results/kpoints/${PREFIX}.out \
         | tail -1 | awk '{print $5}')

  if [ -z "$ETOT" ]; then
    echo "  WARNING: No energy found for k=${K}. Check output."
    ETOT="FAILED"
  fi

  echo "  k-mesh = ${K}x${K}x${K}  ->  Etot = ${ETOT} Ry"
  echo "${K}x${K}x${K}   ${ETOT}" >> results/kpoints/etot_vs_kmesh.dat

done

rm -f tmp_scf.in

echo ""
echo "============================================"
echo " RESULTS saved to: results/kpoints/etot_vs_kmesh.dat"
echo "============================================"
cat results/kpoints/etot_vs_kmesh.dat
echo ""
echo " Choose the k-mesh where |ΔEtot| < 1 meV/atom (~7.4e-5 Ry)"
