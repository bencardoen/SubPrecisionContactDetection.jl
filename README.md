# SubPrecisionContactDetection.jl
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/bencardoen/SubPrecisionContactDetection.jl/tree/main.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/bencardoen/SubPrecisionContactDetection.jl/tree/main) [![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0) [![codecov](https://codecov.io/gh/bencardoen/SubPrecisionContactDetection.jl/branch/main/graph/badge.svg?token=HJ7KIHBZC0)](https://codecov.io/gh/bencardoen/SubPrecisionContactDetection.jl)[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://bencardoen.github.io/SubPrecisionContactDetection.jl/dev/)[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://bencardoen.github.io/SubPrecisionContactDetection.jl/dev/)


Detects sub-precision contacts between subcellular organelles in 2 and 3D STED (precision ~ 50-150nm)
superresolution microscopy, for example endoplasmum reticulum and mitochondria (contacts ~ 20-100nm).

Where a pixel precise segmentation is not feasible due to the precision of the microscope, and colocalization does not describe the interface in a meaningful way, SubPrecisionContactDetection can reconstruct the plausible interface between the organelles.

An example rendering of the postprocessed contact zones (white) between endoplasmum reticulum (green) and mitochondria (red) is shown here [(source)](https://www.biorxiv.org/content/10.1101/2022.06.23.497346v1.full.pdf):

![](example.png)

### Features
- Fast: using multiple threads, and Julia's fast LLVM JIT code
- Reproducible: tests ensure backwards compatibility
- Configurable: Can process deconvolved or raw images, with optional extra denoising
- Rich: provides interpretable features for each detected contact
- Confidence map for each voxel (each voxel has a p-value)


### Cite
If you find this project useful, please cite
```bibtex
@article{cardoen2023membrane,
  title={Membrane contact site detection (MCS-DETECT) reveals dual control of rough mitochondria--ER contacts},
  author={Cardoen, Ben and Vandevoorde, Kurt R and Gao, Guang and Ortiz-Silva, Milene and Alan, Parsa and Liu, William and Tiliakou, Ellie and Vogl, A Wayne and Hamarneh, Ghassan and Nabi, Ivan R},
  journal={Journal of Cell Biology},
  volume={223},
  number={1},
  pages={e202206109},
  year={2023},
  publisher={Rockefeller University Press}
}
```

