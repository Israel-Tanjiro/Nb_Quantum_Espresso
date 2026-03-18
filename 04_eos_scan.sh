#!/bin/bash
# ============================================================
#  EOS SCAN: Equation of State for Nb (BCC)
#  Scans 9 volumes around equilibrium and saves E vs V
#  Run: bash 04_eos_scan.sh
# ============================================================

# --- USER SETTINGS ------------------------------------------
PW_EXEC="pw.x"
PSEUDO_DIR="./pseudo"
TMP_DIR="./tmp"
PSEUDO_FILE="nb_pbe_v1.uspp.F.UPF"
NCORE=6
# ------------------------------------------------------------

# M2 Pro optimized
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p ../../tmp results/eos

# celldm(1) values to scan in Bohr (~±3% around 6.235)
#CELLDMS="6.05 6.10 6.15 6.20 6.235 6.25 6.30 6.35 6.40"
# NEW — uniform 0.05 Bohr spacing, centered on 6.25:
CELLDMS="6.00 6.05 6.10 6.15 6.20 6.25 6.30 6.35 6.40 6.45 6.50"
echo "============================================"
echo " Starting EOS scan for Nb (BCC)"
echo " ecutwfc = 70 Ry | k-mesh = 16x16x16"
echo " Scanning celldm(1) from 6.05 to 6.40 Bohr"
echo "============================================"
echo ""

# Header for results file
echo "# celldm(Bohr)   a(Ang)       V(Ang^3)     Etot(Ry)" \
     > results/eos/eos_data.dat

for CELLDM in $CELLDMS; do

  PREFIX="Nb_eos_${CELLDM}"

  # Convert celldm to Angstrom and Volume for reference
  # a(Ang) = celldm * 0.529177
  # V(Ang^3) = (a/sqrt(2))^3 * 2  for BCC = a^3/2...
  # actually V_BCC = a^3 / 2 but QE celldm(1) is a in Bohr
  # a_ang = celldm * 0.529177
  # V = a_ang^3 / 2  (BCC has 1 atom, primitive cell volume = a^3/2)

  echo "  Running celldm(1) = ${CELLDM} Bohr ..."

  cat > results/eos/${PREFIX}.in << EOF
&CONTROL
  calculation   = 'scf'
  prefix        = '${PREFIX}'
  outdir        = '${TMP_DIR}'
  pseudo_dir    = '${PSEUDO_DIR}'
/
&SYSTEM
  ibrav         = 3
  celldm(1)     = ${CELLDM}
  nat           = 1
  ntyp          = 1
  ecutwfc       = 70
  ecutrho       = 560
  occupations   = 'smearing'
  smearing      = 'marzari-vanderbilt'
  degauss       = 0.03
/
&ELECTRONS
  conv_thr      = 1.0d-10
  mixing_beta   = 0.4
  mixing_mode   = 'plain'
  electron_maxstep = 200
/
ATOMIC_SPECIES
  Nb  92.906  ${PSEUDO_FILE}
ATOMIC_POSITIONS {alat}
  Nb  0.0  0.0  0.0
K_POINTS {automatic}
  16 16 16  0 0 0
EOF

  mpirun -np $NCORE $PW_EXEC -npool $NCORE \
         -in  results/eos/${PREFIX}.in \
         > results/eos/${PREFIX}.out 2>&1

  # Extract total energy
  ETOT=$(grep "!    total energy" results/eos/${PREFIX}.out \
         | tail -1 | awk '{print $5}')

  if [ -z "$ETOT" ]; then
    echo "  WARNING: SCF failed for celldm=${CELLDM} — check output"
    ETOT="FAILED"
    echo "${CELLDM}   FAILED   FAILED   FAILED" >> results/eos/eos_data.dat
  else
    # Compute a in Angstrom and volume in Ang^3
    A_ANG=$(echo "scale=8; $CELLDM * 0.529177" | bc)
    V_ANG=$(echo "scale=8; $A_ANG * $A_ANG * $A_ANG / 2" | bc)
    echo "  celldm = ${CELLDM} Bohr | a = ${A_ANG} Ang | V = ${V_ANG} Ang^3 | E = ${ETOT} Ry"
    echo "${CELLDM}   ${A_ANG}   ${V_ANG}   ${ETOT}" >> results/eos/eos_data.dat
  fi

done

echo ""
echo "============================================"
echo " EOS scan complete!"
echo " Results saved to: results/eos/eos_data.dat"
echo "============================================"
echo ""
cat results/eos/eos_data.dat
echo ""
echo " Next step: python3 fit_eos.py"
