#!/bin/bash
# ================================================================
#  run_overnight.sh
#  Chains the full phonon workflow safely
#  Run before sleep: bash run_overnight.sh
#  Check results in morning: python3 plot_phonon.py
# ================================================================

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

LOG="overnight_run.log"
echo "Started: $(date)" | tee $LOG
echo "" | tee -a $LOG

# Prevent Mac from sleeping during calculation
# (requires caffeinate, built into macOS)
echo "Preventing Mac sleep with caffeinate..."
caffeinate -i &
CAFFEINATE_PID=$!

# Run the full phonon pipeline
bash 05_phonon_pipeline.sh 2>&1 | tee -a $LOG

STATUS=${PIPESTATUS[0]}

# Kill caffeinate when done
kill $CAFFEINATE_PID 2>/dev/null

echo "" | tee -a $LOG
echo "Finished: $(date)" | tee -a $LOG

if [ $STATUS -eq 0 ]; then
    echo "SUCCESS — run python3 plot_phonon.py to see results" | tee -a $LOG
    # Optional: play a sound when done (macOS)
    afplay /System/Library/Sounds/Glass.aiff 2>/dev/null
else
    echo "FAILED with exit code $STATUS — check overnight_run.log" | tee -a $LOG
    afplay /System/Library/Sounds/Basso.aiff 2>/dev/null
fi
