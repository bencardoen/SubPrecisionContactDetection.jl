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
#SBATCH --array=1-CELLS  # change to nr of cells

## Array batch script for use on clusters.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# Copyright 2020-2022, Ben Cardoen

set -euo pipefail

NOW=$(date +"%m_%d_%Y_HH%I_%M")
echo "Starting setup at $NOW"

## Ensure the singularity image is in place
#cp "<location_of_image.sif>" "$SLURM_TMPDIR/mcsdetect.sif"
cp mcsdetect.sif $SLURM_TMPDIR/mcsdetect.sif

echo "Starting task $SLURM_ARRAY_TASK_ID"
IDIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" in.txt)
ODIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" out.txt)

# Test alphas
IMAGE="$SLURM_TMPDIR/mcsdetect.sif"
P="/opt/SubPrecisionContactDetection.jl"


module load singularity
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SINGULARITY_BINDPATH="/scratch/$USER,$SLURM_TMPDIR"
export SINGULARITY_CACHEDIR="$STMP/singularity/cache"

mkdir -p $ODIR
for ALPHA in 0.001 0.01 0.025 0.05
do
    echo "Starting with alpha $ALPHA"
    mkdir -p $ODIR/$ALPHA
    NOW=$(date +"%m_%d_%Y_HH%I_%M")
    singularity exec $IMAGE julia --project=$P --sysimage=$P/sys_img.so $P/scripts/ercontacts.jl  --inpath $IDIR -r "*[1,2].tif" -w 2 --deconvolved --sigmas 2.5-2.5-1.5 --outpath  $ODIR/$ALPHA --alpha $ALPHA --beta  $ALPHA -c 1 -v 2000 --mode=decon 2>&1 | tee -a log_$NOW.txt
done
NOW=$(date +"%m_%d_%Y_HH%I_%M")

echo "DONE at ${NOW}"
