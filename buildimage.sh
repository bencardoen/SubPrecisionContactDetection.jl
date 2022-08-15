#!/bin/bash
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

# Zips the (cloned) repo, builds the Singularity image
echo "This script requires root and internet access to github."
echo "This script assumes /dev/shm is writeable by your account (and exists)"
set -euo pipefail
NOW=$(date +"%m--%d--%Y ~ %I:%M:%S")
echo "Starting processing at $NOW"

CUR=`pwd`
TMP=/dev/shm
cd $TMP
echo "Cloning repo in $TMP"
git clone git@github.com:bencardoen/SubPrecisionContactDetection.jl.git
echo "Creating archive"
zip -rq SubPrecisionContactDetection.jl.zip SubPrecisionContactDetection.jl
rm -rf $TMP/SubPrecisionContactDetection.jl
echo "Cleaning cloned repo from $TMP done"
mv $TMP/SubPrecisionContactDetection.jl.zip $CUR
echo "Archive complete"
cd $CUR

echo "Building Singularity image"
sudo singularity build image.sif ./singularity_recipes/recipe.def
