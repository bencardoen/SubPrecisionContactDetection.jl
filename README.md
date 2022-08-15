# SubPrecisionContactDetection.jl

Detects sub-precision contacts between subcellular organelles in 2 and 3D STED
superresolution microscopy, for example endoplasmum reticulum and mitochondria.

Where a pixel precise segmentation is not feasible due to the precision of the microscope, and colocalization does not describe the interface in a meaningful way, SubPrecisionContactDetection can reconstruct the plausible interace between the organelles.

### Features
- Fast: using multiple threads, and Julia's fast LLVM JIT code
- Reproducible : tests ensure backwards compatibility
- Configurable : Can process deconvolved or raw images, with optional extra denoising
- Rich : provides interpretable features for each detected contact

## Test status
[![CircleCI](https://circleci.com/gh/bencardoen/SubPrecisionContactDetection.jl/tree/main.svg?style=svg&circle-token=8f7bfd5e06262a0eb9003884e1b543dadbbd0e53)](https://circleci.com/gh/bencardoen/SubPrecisionContactDetection.jl/tree/main)

## Code coverage
[![codecov](https://codecov.io/gh/bencardoen/SubPrecisionContactDetection.jl/branch/main/graph/badge.svg?token=V7DB0LTIGI)](https://codecov.io/gh/bencardoen/SubPrecisionContactDetection.jl)

## Installation
[Get Julia](https://julialang.org/).
We currently develop on 1.7-1.8, backwards compatibility may vary.

### Singularity
If you cannot or do not want to install the dependencies yourself, you can use the [Singularity](https://duckduckgo.com/?t=ffab&q=singularity+ce+docs&ia=web) image, which has all dependencies pre-installed.
To run Singularity on Windows, set up [WSL2](https://www.blopig.com/blog/2021/09/using-singularity-on-windows-with-wsl2/).
```bash
singularity pull library://bcvcsert/mcsdetect/mcsdetect_f35_j1.7:0.0.3
# OR
singularity pull library://bcvcsert/mcsdetect/mcsdetect_f35_j1.6:0.0.3
# then
singularity exec mcsdetect.sif julia --project=/opt/SubPrecisionContactDetection.jl -e 'your code'
# or
singularity shell mcsdetect.sif julia --project=/opt/SubPrecisionContactDetection.jl # Interactive
```
#### Optimized version
You can use the pre-compiled version to get a significant boost in execution speed:
```
singularity exec mcsdetect.sif julia --project=/opt/SubPrecisionContactDetection.jl --sysimage=/opt/SubPrecisionContactDetection.jl/sys_img.so -e 'using SubPrecisionContactDetection; SubPrecisionContactDetection.spear(zeros(1024, 1024), zeros(1024, 1024))'
```

### Install as a global package
```julia
julia
julia> using Pkg;
julia> Pkg.add(url="https://github.com/bencardoen/ERGO.jl.git")
julia> Pkg.add(url="https://github.com/bencardoen/SPECHT.jl.git")
julia> Pkg.add(url="https://github.com/bencardoen/SubPrecisionContactDetection.jl.git")
julia> Pkg.activate(".")
julia> Pkg.instantiate(".")
julia> Pkg.build(".")
julia> Pkg.test("SubPrecisionContactDetection")
```
ERGO and SPECHT are not yet part of the Julia package collection, so they won't be automatically pulled in, which is why we're first adding them.

** Note ** SubPrecisionContactDetection relies on 2 Python packages, which we install for you (with a separate Python installation). If you prefer to use your own venv/Conda, remove the ENV["PYTHON"] step (not recommended unless you know what you're doing).

### Install locally (to use the processing scripts)
```bash
git clone https://github.com/bencardoen/SubPrecisionContactDetection.jl.git
cd SubPrecisionContactDetection.jl
```
```julia
using Pkg; Pkg.activate("."); Pkg.build(); Pkg.instantiate(); Pkg.test();
```

## Detect contacts
```bash
julia --project=. ./src/ercontacts.jl --inpath ./in -r "*[1,2].tif" -w 2 --sigmas 2.5-2.5-1.5 --outpath  ./out --alpha 0.01 --beta 0.01 -c 1 -v 2000 --mode=non-decon 2>&1 | tee -a log_test.txt
```
Where

* --project=. : ensures you run in the local cloned repo
* --{in|out}path : directories where tif files can be found
* -r : regex to tif files, e.g. *[1,2].tif indicates channel 1 and 2 will have filenames ending in 1,2.tif respectively.
* -w : windowsize, >1 or higher
* --sigmas : smoothing Gaussian, set < precision
* --alpha : max false positive rate (p-value)
* --beta : max false negative rate (stat. power)
* -c 1: postprocess channel 1
* -v 2000: drop all contacts touching objects in channel 1 with volume < 2000
* --mode=non-decon : input are non deconvolved tiff files
* 2>&1 | tee -a log_test.txt : save any output to log.txt (in addition to showing it in stdout)

Input is expected to be 2x 3D volumes (in tif format), of equal dimensions.

!Do not use deconvolved images unless you really know what you're doing.!

#### Output
- skeleton_contacts.tif
- channel_[1,2].tif
- (non)_vesicle_contacts: mitochondria (channel 1) with volume < 2000 are considered vesicles, split contacts so you can visualize them separately
- channel_1_(non_)vesicle.tif : channel 1 objects (mitochondria) split into < 2000 and > 2000 objects
- raw|gradient|eroded.tif : stages of progressively computed contacts, all but 'eroded' are debug output
- *.csv : features for each contact

### Notes
- RAM usage : 500x500x20 : ~ 5GB RAM, 2000x2000x70: ~ 50GB RAM, and so on.
- Will use all of JULIA_NUM_THREADS to run in parallel. > 8 is overkill, so set to 4-8 at most: ```export JULIA_NUM_THREADS=4```

### Cite
If you find this project useful, please cite
```bibtex
```

### Cluster scripts
In folder hpcscripts, you'll find 3 scripts, intended to make cluster processing easier.
* sbatch.sh : executes with preset parameters the detection on 1 cell
* array_sbatch.sh : executes any number of cells in parallel
* buildfilelist.sh : builds the list of cells to process for array_sbatch.sh

Usage
```bash
./buildfilelist.sh indir outdir
# 2 files, inlist.txt and outlist.txt are generated
```

Get the number of cells
```bash
wc inlist.txt| awk '{print $1}'
```
In the array_sbatch script, change line
```bash
#SBATCH --array=1-<NROFCELLS>
```
NOTE: Change the account and email fields!!!

Submit to SLURM Scheduler
```
sbatch array_sbatch.sh
```
That's it.
