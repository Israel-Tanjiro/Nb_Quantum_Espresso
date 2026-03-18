#!/bin/bash
# ================================================================
#  PHONON PIPELINE FOR Nb (BCC) — Complete Workflow
#  SCF → ph.x → q2r → matdyn → plot
#
#  Converged parameters:
#    Pseudopotential : nb_pbe_v1.uspp.F.UPF (GBRV, US, PBE, 13e)
#    ecutwfc         : 70 Ry
#    ecutrho         : 560 Ry
#    k-mesh          : 16x16x16
#    degauss         : 0.03 Ry
#    celldm(1)       : 6.2521 Bohr  (DFT optimized)
#
#  Run: bash 05_phonon_pipeline.sh
#  Estimated time on M2 Pro 32GB (6 cores):
#    SCF    :  ~10 min
#    ph.x   :  ~2-4 hrs  (4x4x4 q-mesh)
#    q2r    :  ~1 min
#    matdyn :  ~1 min
#    plot   :  ~1 min
# ================================================================

# ---- USER SETTINGS -------------------------------------------
PW_EXEC="pw.x"
PH_EXEC="ph.x"
Q2R_EXEC="q2r.x"
MATDYN_EXEC="matdyn.x"
NCORE=6
PSEUDO_DIR="./pseudo"
TMP_DIR="./tmp"
PSEUDO_FILE="nb_pbe_v1.uspp.F.UPF"
RESULTS_DIR="results/phonon"
# --------------------------------------------------------------

# M2 Pro optimized thread settings
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p $TMP_DIR $RESULTS_DIR

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo ""
echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}   PHONON DISPERSION PIPELINE FOR Nb (BCC)${NC}"
echo -e "${BLUE}================================================================${NC}"
echo -e "   Cores    : ${NCORE}"
echo -e "   ecutwfc  : 70 Ry"
echo -e "   k-mesh   : 16x16x16"
echo -e "   q-mesh   : 4x4x4"
echo -e "   celldm(1): 6.2521 Bohr"
echo -e "${BLUE}================================================================${NC}"
echo ""

# ================================================================
#  STEP 1 — SCF calculation
# ================================================================
echo -e "${YELLOW}[STEP 1/4] SCF calculation ...${NC}"
START_SCF=$SECONDS

cat > $RESULTS_DIR/Nb_scf.in << EOF
&CONTROL
  calculation   = 'scf'
  prefix        = 'Nb'
  outdir        = '${TMP_DIR}'
  pseudo_dir    = '${PSEUDO_DIR}'
/
&SYSTEM
  ibrav         = 3
  celldm(1)     = 6.2521
  nat           = 1
  ntyp          = 1
  ecutwfc       = 70.0
  ecutrho       = 560.0
  occupations   = 'smearing'
  smearing      = 'marzari-vanderbilt'
  degauss       = 0.03
/
&ELECTRONS
  conv_thr      = 1.0d-12
  mixing_beta   = 0.4
  mixing_mode   = 'plain'
  electron_maxstep = 200
/
ATOMIC_SPECIES
  Nb  92.906  ${PSEUDO_FILE}
ATOMIC_POSITIONS {alat}
  Nb  0.0  0.0  0.0
K_POINTS {automatic}
  12 12 12  0 0 0
EOF

mpirun -np $NCORE $PW_EXEC -npool $NCORE \
       -in  $RESULTS_DIR/Nb_scf.in \
       > $RESULTS_DIR/Nb_scf.out 2>&1

# Check SCF success
if grep -q "convergence has been achieved" $RESULTS_DIR/Nb_scf.out; then
    ETOT=$(grep "!    total energy" $RESULTS_DIR/Nb_scf.out | tail -1 | awk '{print $5}')
    SCF_TIME=$(( SECONDS - START_SCF ))
    echo -e "   ${GREEN}SCF converged!${NC}  Etot = ${ETOT} Ry  (${SCF_TIME}s)"
else
    echo -e "   ${RED}ERROR: SCF did not converge — check $RESULTS_DIR/Nb_scf.out${NC}"
    exit 1
fi

echo ""
# 
# # ================================================================
# #  STEP 2 — Phonon calculation (ph.x)
# # ================================================================
# echo -e "${YELLOW}[STEP 2/4] Phonon calculation (ph.x) — this takes ~2-4 hrs ...${NC}"
# echo -e "   Started at: $(date)"
# echo -e "   You can monitor progress with:"
# echo -e "   tail -f $RESULTS_DIR/Nb_ph.out"
# echo ""
# START_PH=$SECONDS
#
# cat > $RESULTS_DIR/Nb_ph.in << EOF
# Phonon calculation for Nb BCC
#  &INPUTPH
#   prefix       = 'Nb'
#   outdir       = '${TMP_DIR}'
#   fildyn       = '${RESULTS_DIR}/Nb.dyn'
#   ldisp        = .true.
#   nq1          = 4
#   nq2          = 4
#   nq3          = 4
#   tr2_ph       = 1.0d-14
#   alpha_mix(1) = 0.3
#   recover      = .true.
# /
# EOF
#
# mpirun -np $NCORE $PH_EXEC -npool $NCORE \
#        -in  $RESULTS_DIR/Nb_ph.in \
#        > $RESULTS_DIR/Nb_ph.out 2>&1
#
# # Check ph.x success
# if grep -q "PHonon code" $RESULTS_DIR/Nb_ph.out && \
#    ! grep -q "stopping" $RESULTS_DIR/Nb_ph.out; then
#     PH_TIME=$(( SECONDS - START_PH ))
#     PH_MIN=$(( PH_TIME / 60 ))
#     echo -e "   ${GREEN}Phonon calculation complete!${NC}  (${PH_MIN} min)"
# else
#     # Check for normal completion differently
#     if grep -q "PHONON" $RESULTS_DIR/Nb_ph.out; then
#         PH_TIME=$(( SECONDS - START_PH ))
#         PH_MIN=$(( PH_TIME / 60 ))
#         echo -e "   ${GREEN}Phonon calculation complete!${NC}  (${PH_MIN} min)"
#     else
#         echo -e "   ${RED}ERROR: ph.x may have failed — check $RESULTS_DIR/Nb_ph.out${NC}"
#         echo -e "   ${YELLOW}Attempting to continue with q2r anyway ...${NC}"
#     fi
# fi
#
# echo ""
#
# # ================================================================
# #  STEP 3 — q2r: dynamical matrices → interatomic force constants
# # ================================================================
# echo -e "${YELLOW}[STEP 3/4] q2r — Fourier transform to real space IFCs ...${NC}"
#
# cat > $RESULTS_DIR/Nb_q2r.in << EOF
#  &INPUT
#   fildyn  = '${RESULTS_DIR}/Nb.dyn'
#   zasr    = 'simple'
#   flfrc   = '${RESULTS_DIR}/Nb.fc'
#  /
# EOF
#
# $Q2R_EXEC < $RESULTS_DIR/Nb_q2r.in > $RESULTS_DIR/Nb_q2r.out 2>&1
#
# if grep -q "q2r" $RESULTS_DIR/Nb_q2r.out || [ -f "${RESULTS_DIR}/Nb.fc" ]; then
#     echo -e "   ${GREEN}q2r complete!${NC}  IFC file: $RESULTS_DIR/Nb.fc"
# else
#     echo -e "   ${RED}ERROR: q2r failed — check $RESULTS_DIR/Nb_q2r.out${NC}"
#     exit 1
# fi
#
# echo ""
#
# # ================================================================
# #  STEP 4 — matdyn: phonon dispersion along high-symmetry path
# #  BCC high-symmetry path: Γ → H → P → Γ → N → Σ → Γ
# # ================================================================
# echo -e "${YELLOW}[STEP 4/4] matdyn — phonon dispersion along BCC path ...${NC}"
#
# cat > $RESULTS_DIR/Nb_matdyn.in << EOF
#  &INPUT
#   asr      = 'simple'
#   flfrc    = '${RESULTS_DIR}/Nb.fc'
#   flfrq    = '${RESULTS_DIR}/Nb.freq'
#   flvec    = '${RESULTS_DIR}/Nb.modes'
#   fleig    = '${RESULTS_DIR}/Nb.eig'
#   q_in_band_form = .true.
#   dos      = .false.
#  /
# 8
#   0.000  0.000  0.000  30   ! Gamma
#   1.000  0.000  0.000  25   ! H
#   0.500  0.500  0.500  25   ! P
#   0.000  0.000  0.000  30   ! Gamma
#   0.000  1.000  0.000  25   ! N
#   0.500  0.500  0.000  20   ! Sigma (midpoint)
#   0.000  0.000  0.000  10   ! Gamma (again)
#   0.500  0.500  0.000   1   ! end
# EOF
#
# $MATDYN_EXEC < $RESULTS_DIR/Nb_matdyn.in > $RESULTS_DIR/Nb_matdyn.out 2>&1
#
# if [ -f "$RESULTS_DIR/Nb.freq" ]; then
#     echo -e "   ${GREEN}matdyn complete!${NC}  Frequencies: $RESULTS_DIR/Nb.freq"
# else
#     echo -e "   ${RED}ERROR: matdyn failed — check $RESULTS_DIR/Nb_matdyn.out${NC}"
#     exit 1
# fi
#
# # ================================================================
# #  STEP 5 — Generate band path labels using bands_x helper
# # ================================================================
# echo ""
# echo -e "${YELLOW}[STEP 5/5] Preparing dispersion data for plotting ...${NC}"
#
# # Use plotband.x to generate gnuplot-friendly file
# cat > $RESULTS_DIR/Nb_plotband.in << EOF
# ${RESULTS_DIR}/Nb.freq
# 0 500
# ${RESULTS_DIR}/Nb_bands.plot
# ${RESULTS_DIR}/Nb_bands.ps
# 0.0
# 60.0 0.0
# EOF
#
# # Also run matdyn for phonon DOS
# cat > $RESULTS_DIR/Nb_dos.in << EOF
#  &INPUT
#   asr      = 'simple'
#   flfrc    = '${RESULTS_DIR}/Nb.fc'
#   flfrq    = '${RESULTS_DIR}/Nb_dos.freq'
#   dos      = .true.
#   fldos    = '${RESULTS_DIR}/Nb.dos'
#   nk1      = 20
#   nk2      = 20
#   nk3      = 20
#   deltaE   = 1.0
#  /
# EOF
#
# $MATDYN_EXEC < $RESULTS_DIR/Nb_dos.in > $RESULTS_DIR/Nb_dos.out 2>&1
#
# if [ -f "$RESULTS_DIR/Nb.dos" ]; then
#     echo -e "   ${GREEN}Phonon DOS complete!${NC}"
# fi
#
# # ================================================================
# #  Summary
# # ================================================================
# TOTAL_TIME=$(( SECONDS ))
# TOTAL_MIN=$(( TOTAL_TIME / 60 ))
#
# echo ""
# echo -e "${BLUE}================================================================${NC}"
# echo -e "${GREEN}   PHONON PIPELINE COMPLETE!${NC}"
# echo -e "${BLUE}================================================================${NC}"
# echo -e "   Total time     : ${TOTAL_MIN} min"
# echo -e "   Dispersion data: $RESULTS_DIR/Nb.freq"
# echo -e "   Phonon DOS     : $RESULTS_DIR/Nb.dos"
# echo -e "   Log files      : $RESULTS_DIR/*.out"
# echo ""
# echo -e "   Next step — plot results:"
# echo -e "   ${YELLOW}python3 plot_phonon.py${NC}"
# echo -e "${BLUE}================================================================${NC}"
# echo ""
