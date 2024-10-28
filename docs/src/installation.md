# Installation
This project is developed using [Julia](https://julialang.org/).

For ease of use and to maximize reproducibility we also provide container images using Singularity.

## Source code
The fastest and easiest way is to clone the repository using [Git](https://git-scm.com/). Alternatively you can download zipped [releases](https://github.com/bencardoen/SubPrecisionContactDetection.jl/archive/refs/heads/main.zip).
The below assumes you have Julia 1.9 or higher installed.

In an empty, new folder:
```bash
git clone https://github.com/bcardoen/SubPrecisionContactDetection.jl.git
cd SubPrecisionContactDetection.jl
julia --project=. -e `using Pkg; Pkg.build(); Pkg.test()`
```
This downloads the source code in a subfolder, builds it with all dependencies and tests it.

Once you have this, you can use either a terminal or an IDE (e.g. Visual Studio Code) to work with the source code to process new datasets.

## Singularity/Apptainer
We provide a preconfigured container image with all dependencies [here](https://cloud.sylabs.io/library/bcvcsert/subprecisioncontactdetection/mcsdetect).

