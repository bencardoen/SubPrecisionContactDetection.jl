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

### Updating
If you want to get the latest version, using git do:
```bash
git pull origin main
```
Then make sure everything still works
```
julia --project=. -e 'using Pkg; Pkg.build();Pkg.test();`
```

This assumes you did not switch brances or modified the code. 

## Singularity/Apptainer
We provide a preconfigured container image with all dependencies [here](https://cloud.sylabs.io/library/bcvcsert/subprecisioncontactdetection/mcsdetect).


## Dependencies
All dependencies are automatically installed, however, you may run into issues if you have a non-standard Python installation. 
See the [build](https://github.com/bencardoen/SubPrecisionContactDetection.jl/build/build.jl) script for details, but when in doubt, set the environment variable PYTHON to either your Python installation or the empty string and rebuild
```bash
export PYTHON=""
julia --project=. -e 'using Pkg; Pkg.build'
```

!!! note "Attention"
    For the remainder of this document we assume all commands are run inside the cloned directory, e.g. `SubPrecisionContactDetection.jl`.