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
# Copyright 2021-2022, Ben Cardoen

export LOCALPKG=/opt/DataCurator.jl
export JLMJV=1.8
export JLV=$JLMJV.5
export PATH=/opt/julia/julia-$JLV/bin:$PATH
export JULIA_DEPOT_PATH=/opt/juliadepot

echo "Starting ..."

julia --project=/opt/SubPrecisionContactDetection.jl --sysimage=/opt/SubPrecisionContactDetection.jl/sys_img.so "$@"
