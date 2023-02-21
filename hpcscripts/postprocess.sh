# Copyright (C) 2023 Ben Cardoen
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!/bin/bash
#SBATCH --account=ACCOUNT
#SBATCH --mem=120G
#SBATCH --cpus-per-task=6
#SBATCH --time=18:00:00
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

set -euo pipefail

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

NOW=$(date +"%m_%d_%Y_HH%I_%M")
echo "Starting setup at $NOW"

## Ensure the singularity image is in place
#cp "<location_of_image.sif>" "$SLURM_TMPDIR/mcsdetect.sif"
cp mcsdetect.sif $SLURM_TMPDIR/mcsdetect.sif

IDIR=INPUT
ODIR=POSTOUTPUT

# Test alphas
IMAGE="$SLURM_TMPDIR/mcsdetect.sif"
LSRC="/opt/SubPrecisionContactDetection.jl"

module load singularity
export SINGULARITY_BINDPATH="/scratch/$USER,$SLURM_TMPDIR"
singularity exec $IMAGE julia --project=$LSRC --sysimage=$LSRC/sys_img.so $LSRC/scripts/run_cube_sampling_on_dataset.jl  --inpath $IDIR --outpath  $ODIR 2>&1 | tee -a log_$NOW.txt

NOW=$(date +"%m_%d_%Y_HH%I_%M")

echo "DONE at ${NOW}"
