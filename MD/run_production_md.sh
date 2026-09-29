#!/bin/bash
# =============================================================
# Rosavin / rosin / L-arabinose + H2O2 -- production MD
# Continues from a pre-equilibrated NPT checkpoint.
# =============================================================
# Requirements: GROMACS
#   (On Kaggle, GROMACS is not preinstalled; install via conda-forge/bioconda,
#    e.g. with condacolab + `mamba install -c conda-forge -c bioconda gromacs`)
#
# Usage:
#   ./run_production_md.sh <compound> <rep>
#   e.g. ./run_production_md.sh l_arabinose 2
#
# Expects, in md/<compound>/rep<rep>/:
#   md.mdp            production MD parameters
#   npt.gro           structure after NPT equilibration
#   npt.cpt           NPT checkpoint (continuation state)
#   topol.top         topology (includes solute + H2O2 + water/ion itp files)
#   *.itp             solute, H2O2, and position-restraint include files
#                      referenced by topol.top
# =============================================================

set -e  # stop on first error

COMPOUND=$1   # rosavin | rosin | l_arabinose
REP=$2        # 1 | 2

if [ -z "$COMPOUND" ] || [ -z "$REP" ]; then
  echo "Usage: ./run_production_md.sh <compound> <rep>"
  exit 1
fi

RUN_DIR="md/${COMPOUND}/rep${REP}"
cd "$RUN_DIR"

echo "Running production MD for ${COMPOUND}, replicate ${REP} in ${RUN_DIR}"

# ---- 1. Build the production-run input (.tpr), continuing from NPT ------
gmx grompp \
  -f md.mdp \
  -c npt.gro \
  -t npt.cpt \
  -p topol.top \
  -o md.tpr \
  -maxwarn 2

# ---- 2. Run production MD ------------------------------------------------
# -cpi md.cpt lets this be re-run safely to resume from the last checkpoint
# (e.g. after a Kaggle/Colab session timeout) rather than starting over.
gmx mdrun -deffnm md -v -cpi md.cpt

echo "Done. Production run output: ${RUN_DIR}/md.gro, md.xtc, md.log, md.edr"
