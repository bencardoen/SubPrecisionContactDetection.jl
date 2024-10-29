# Postprocessing

Once the contact maps have been computed, you often need quantification and additional filtering. 
For example, coverage, features descriptors, and so forth.

There are three key processing steps disjoint from the actual algorithm output:
- bleedthrough filter
- CSV curation
- Sampling

In the remaining of this document, let us assume `DIR` is the directory where the algorithm saved its output on a full dataset.

## Bleedthrough filter
The background filter removes ghost effects (bleedthrough).
It is run as part of the pipeline, but you can invoke it separately.

!!! note "This is Optional"
    This is entirely optional, but useful if you want to optimize this filter independently.

Suppose we want to filter all tif files ending with "1.tif" or "2.tif" , for z=1 to 1.1 in 0.25 steps, and then compute the object properties.
```julia
julia --project=.  scripts/segment.jl --inpath mydir -z 1.0 -Z 1.1 -s 0.25 -r "*[1,2].tif"
```

For each file, for each filtered value, it will generate:
- `mask.tif`
- `masked.tif`

For all the files, it will generate a CSV with columns, where each row is an object:

- `size` (nr of non zero voxels)
- `weighted` (intensity sum of objects)
- `minimum, Q1, mean, mediam, Q3, maximum, std, kurtosis` : describes the intensity distribution of the object
- `xyspan` : major axis of 2D projection, pixels
- `zrange` : extent in Z
- `zmidpoint` : midpoint in Z slice
- `distance_to_centroid`: distance of this object's centroid to centroid of all objects, describes clustering
- `distance_to_centroid_normalized`: rescale the above to 0-1, where 0 is centroid of all objects, 1 is maximum dispersion
- `centroid_channel_{x,y,z}` : in pixel coordinates, records the centroid of _all_ objects, this is the reference for the distance computation
- `centroid_object_{x,y,z}` : in pixel coordinates, records the object centroid
- `filename` : the source tif file name
- `z` : the z value used
- `eig1-3`: PCA eigenvalues, the can be used for shape descriptors.
- `eig1-3normalized`: eigenvalues rescaled to 1.


!!! warning "Shape features"
    The ``\lambda`` values are disabled by default due given that for very large objects they can stall the pipeline (1e6 voxels).

## CSV Curation
You can run our Python script to aggregate and curate the processed CSV files.

```python
python3 scripts/csvcuration.py --inputdirectory <where you saved the output> --outputdirectory <where you want the new CSV files saved>
```
By default this will look for output produced with ``\alpha`` 0.05, you can override this as needed with `--alpha 0.01` for example.

This will produce:

```
contacts_aggregated.csv             # Contacts aggregated per cell, so 1 row = 1 cell, use this for e.g. mean height, Q95 Volume
contacts_filtered_novesicles.csv    # All contacts, without vesicles
contacts_unfiltered.csv             # All contacts, no filtering
```

## Sampling contacts
In [scripts/run_cube_sampling_on_dataset.jl](https://github.com/bencardoen/SubPrecisionContactDetection.jl/scripts/run_cube_sampling_on_dataset.jl) you'll find a script that samples contacts with a sliding window, to avoid long tail statistics dominating the conclusion of any analysis. The paper goes into more depth why this is beneficial.

```julia
julia --project=. scripts/run_cube_sampling_on_dataset.jl  --inpath DIR --outpath  <where to save your output>
```

A convenience script is provided to further aggregate the output of this stage.

```python
python3 scripts/coverage.py  --inputdirectory DIR --outputdirectory <where to save your ouput>
```

This will print summary output and save a file `coverage_aggregated.csv`. The columns Coverage % mito by contacts, mean per cell and ncontacts mean are the columns you'll be most interested in.

They report the coverage of contacts on mitochondria (minus MDVs), and the number of contacts per sliding window of 5x5x5 voxels.
    
!!! warning "Singularity/Apptainer usage"
    See the [tutorial](https://github.com/NanoscopyAI/tutorial_mcs_detect) on how the run these with a container image. The code is the same, you just need to define some extra variables. 