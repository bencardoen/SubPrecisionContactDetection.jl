# Postprocessing

Once the contact maps have been computed, you often need quantification and additional filtering. 
For example, coverage, features descriptors, and so forth.

## Aggregating CSV files


## Sampling contacts
In [scripts/run_cube_sampling_on_dataset.jl](https://github.com/bencardoen/SubPrecisionContactDetection.jl/scripts/run_cube_sampling_on_dataset.jl) you'll find a script that samples contacts with a sliding window, to avoid long tail statistics dominating the conclusion of any analysis. The paper goes into more depth why this is beneficial.



## Preprocessing and filtering
The background filter removes ghost effects (bleedthrough).
If you want to tune this without invoking the full pipeline, you can do so:

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
    
