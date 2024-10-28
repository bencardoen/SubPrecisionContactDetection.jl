### Output
A brief overview of the generated output follows here:

**Files are saved with the prefix of the original**: If you have files, say, ABC1.tif, and ABC2.tif, then output will be saved ABC_xyz.tif, and so forth.
- $PREFIX_confidence_map.tif: Each voxel holds the confidence (significance, p-value [0-1]) of the corresponding contact.
- $PREFIX__channel_[1,2].tif: The filtered mitochondria and ER channel. Use this to inspect if the z-filter removed too much or too little
- $PREFIX__pre_split_raw.tif: Unfiltered contacts. At voxel x, y, z this will have a correlation value, and its significance can be found in confidence_map[x,y,z].
- $PREFIX_pre_split_gradient.tif: Contacts with the gradient filter applied.
- $PREFIX_split_eroded.tif: Erodes singular voxels that are below the precision of the system
- $PREFIX_3_eroded_volumes_nonsplit.csv : Features computed on the contacts
- $PREFIX_C1_objects.csv: A CSV file where each row describe the features of the segmented objects of that channel (1). So if the code ran on 01.tif and 02.tif, C1 will map to 01.tif, filtered, then processed.
- $PREFIX_C2_objects.csv: A CSV file where each row describe the features of the segmented objects of that channel (1). So if the code ran on 01.tif and 02.tif, C1 will map to 01.tif, filtered, then processed.
See also the postprocessing for further output.
The remainder are debugging outputs that can be traced to their source code in the script.
