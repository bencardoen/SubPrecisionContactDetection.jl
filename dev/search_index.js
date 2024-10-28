var documenterSearchIndex = {"docs":
[{"location":"parameters/#Parameter-selection-and-tuning.","page":"Parameter selection and tuning","title":"Parameter selection and tuning.","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"This section describes the usage of parameters for the scripts and the algorithm.","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Pages = [\"parameters.md\"]\nDepth = 5","category":"page"},{"location":"parameters/#Checking-default-values","page":"Parameter selection and tuning","title":"Checking default values","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"You can easily check the default values and the parameter names by querying the code.","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"params = get_defaults()\nfor k in keys(params)\n    @info \"$k ==> $(params[k])\"\nend","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"This would give:","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"[ Info: lpsigmas → 1-1-1\n[ Info: inregex → *[1,2].tif\n[ Info: minzslice → 1\n[ Info: mode → non-decon\n[ Info: zscore → 3\n[ Info: outpath → \n[ Info: radius → false\n[ Info: sphericity → 1.0\n[ Info: volumethreshold → 0\n[ Info: cube-vesicle-intensity-mean → 0.2\n[ Info: nooutput → false\n[ Info: normalize → false\n[ Info: alpha → 0.05\n[ Info: prc → 1.0\n[ Info: cube-vesicle-sample-size → 5\n[ Info: denominator → 1.0\n[ Info: dimension → 3\n[ Info: weighted → false\n[ Info: beta → 0.05\n[ Info: volumethresholdchannel → 1\n[ Info: noutput → false\n[ Info: save-numerical-data → false\n[ Info: windowsize → 1\n[ Info: skipstats → false\n[ Info: sigmas → 1-1-1\n[ Info: inpath → \n[ Info: filtermode → arithmetic\n[ Info: cube-vesicle-size-ln → 9\n[ Info: dry-run → false\n[ Info: deconvolved → true","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Critical ones are :","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"inpath: Directory with your data\noutpath: Where output should be written\ninregex: \"*[1,2].tif\" Looks for files ending with 1 or 2 (extension tif). If you want to use 4 channels out of 8, and only odd ones, you could write \"*[0,2,4,6].tif\". The pipeline will automatically combine them for you. frac4(2)=6 combinations would be generated, in order, e.g. \"0–2\", \"0–4\", etc.\ndimension : 2 or 3\nmode : non-decon, decon, or both. Do not use non-deconvolved data, stick to 'decon'\nalpha/beta: Significance and power used in probabilistic filtering\nwindowsize: w in k-D means a window of (1 + 2 times w)^k voxels. \nzscore: see below","category":"page"},{"location":"parameters/#Parameter-selection","page":"Parameter selection and tuning","title":"Parameter selection","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"MCS-Detect has multiple parameters that will determine the precision and recall of the predicted contacts.  While a full discussion is available in the paper, here we will give a brief explanation and guidance as to how to set them.","category":"page"},{"location":"parameters/#Z-filter-(background-removal)","page":"Parameter selection and tuning","title":"Z-filter (background removal)","text":"","category":"section"},{"location":"parameters/#Concept","page":"Parameter selection and tuning","title":"Concept","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Because 3D STED has anisotropic resolution in Z (worse in Z than in X/Y), it is possible to see intensity bleedthrough or shadowing across Z.  For example, say you have a mitochondrial vesicle at Z-slice 5.  Bleedthrough can lead to intensity mimicking a faint object at Z-slice 8. The Z-filter removes this by filtering the intensity distribution, per channel. If you set Z=1, all intensity below mu + 1 * sigma is set to zero.","category":"page"},{"location":"parameters/#Guidance","page":"Parameter selection and tuning","title":"Guidance","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"A z-value is that is too high will cause false negatives because you're removing intensity from the organelles, not the background. A too low value will included possible contacts between organelles and phantom intensity, e.g. false positives. A value of z=3 is used for the paper, derived from the size of the cell and the anisotropy.  Recommended usage is to test Z-values on a single representative cell, and plot the organelle volume, in combination with visual inspection.  Instructions on how to do this and accompanying scripts can be found here.","category":"page"},{"location":"parameters/#Window-size-(w)","page":"Parameter selection and tuning","title":"Window size (w)","text":"","category":"section"},{"location":"parameters/#Concept-2","page":"Parameter selection and tuning","title":"Concept","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Correlation requires a comparison between two vectors of data, in 2- or 3D images this means a window size.  If you set w=2 the window will be (2*2+1)^D for D dimensions. So 5x5 in 2D, 5x5x5 in 3D. w=1 would be 3x3, or 3x3x3 and so forth.","category":"page"},{"location":"parameters/#Guidance-2","page":"Parameter selection and tuning","title":"Guidance","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"A too large value will consume more memory, and will miss finer patterns.  A too small value will fail to capture large patterns.  So what, then is 'too' small or large? At a minimum, the window should cover the width of the contact, but no more than 2.5x.  The interested reader will detect similarities with how resolution and pixel-dimensions relate. I will give an example to give a more actionable insight: Let us assume pixel precision is 50nm in X, Y, and 75nm in Z. Say the expected contacts you wish to capture are 0-25nm.  In this case w=1 would be sufficient, because a window of 3x3x3 would span 150nm lateral, and 225nm axial.  W=2 would mean 250nm lateral and 375nm axial, which is likely too large, it would be dominated by differentials that are unrelated to the contact. Important The window size determines the statistical power of the correlation. A 3x3 window in 2D has limited statistical power. See below.","category":"page"},{"location":"parameters/#Alpha-and-Beta,-or-confidence,-significance,-and-power.","page":"Parameter selection and tuning","title":"Alpha and Beta, or confidence, significance, and power.","text":"","category":"section"},{"location":"parameters/#Concept-3","page":"Parameter selection and tuning","title":"Concept","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"A correlation is a statistical estimator, and comes with a confidence value ('p-value').  Alpha control what acceptable levels of confidence are allowed, whereas beta controls statistical power.  A recap from statistics:","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Significance (alpha): The probability that an observed difference is not due to random effects\nPower (beta): The probability that you can observe a given difference (of a given magnitude)","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"What does this mean in practice? We can compute what is minimal observable correlation you can detect, given alpha and beta. First, the 2D case (so 3x3, 5x5, ...)","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"(Image: minr2d.png)","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Next, 3D:","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"(Image: minr3d.png)","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Trouble reading these plots? Let's say you use a 3x3x3 window (w=1, in 3D).  If you set alpha=beta=0.05 (95% confidence and power), then the smallest possible observable correlation is 0.665. (In the 2nd plot, X=27, Y=0.665).","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Suppose you increase the window to w=2, 3D, then you have 0.341 (X=125, Y=0.341).","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"If you want to have the same minimum correlation in 3D with a window of 27, you would need to change your alpha and beta to 0.35","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"We can also plot this ","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"(Image: minrkd.png)","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"The functions to compute this are available for you as well:","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"# w=2, 2D\nminr = compute_min_r_for_sample_corr(25, 0.05, 0.05)","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"and","category":"page"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"# r=0.2, 2D\nwindow = compute_sample_size_for_min_corr(0.2, 0.05, 0.05)","category":"page"},{"location":"parameters/#Guidance-3","page":"Parameter selection and tuning","title":"Guidance","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"If you keep the window the same, and go to 2D, set alpha and beta from to have the same recall.\nIf precision is too low, reduce alpha and beta (e.g. 0.05 to 0.1, or 0.25).\nIf recall is too high (artifacts), increase alpha and beta (0.05 to 0.01 or 0.001)","category":"page"},{"location":"parameters/#Vesicle-filtering","page":"Parameter selection and tuning","title":"Vesicle filtering","text":"","category":"section"},{"location":"parameters/#Concept-4","page":"Parameter selection and tuning","title":"Concept","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"The postprocessing scripts use size (logarithm) and mean intensity of vesicles to filter them out.  This can only be empirically estimated. ","category":"page"},{"location":"parameters/#Guidance-4","page":"Parameter selection and tuning","title":"Guidance","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Plot the intensity and object sizes of the mitochondria channel, and look for a separation between large and bright objects, versus small and faint. Off the shelf clustering methods can be of help. Alternatively, segment the image before processing.  NOTE Contact detection does not differentiate between mitochondria and vesicles, the interaction may be functionally different, but the contacts are no less real.","category":"page"},{"location":"parameters/#Sampling","page":"Parameter selection and tuning","title":"Sampling","text":"","category":"section"},{"location":"parameters/#Concept-5","page":"Parameter selection and tuning","title":"Concept","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"Because contacts are large and infrequent, or small and frequent, the statistical analysis can be unstable.  More precisely, the distribution is long tailed containing extreme values, and those extreme values are often the ones of interest (e.g. ribomerc). To offset this, the coverage computation and local density (nr of contacts/window) uses a window, defaulting to 5x5x5 (which corresponds to w=2). ","category":"page"},{"location":"parameters/#Guidance-5","page":"Parameter selection and tuning","title":"Guidance","text":"","category":"section"},{"location":"parameters/","page":"Parameter selection and tuning","title":"Parameter selection and tuning","text":"The smaller you set this, the more you split objects apart.  Ideally you set this window to be no smaller than the largest expected object. Sampling windows do not overlap, and mitochondria that are only partially visible in a window (few voxels), are discarded.","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The below assumes you have the source code installed as per the installation instructions.","category":"page"},{"location":"tutorial/#Processing-a-single-cell-(3D-STED)-with-two-channels","page":"Tutorial","title":"Processing a single cell (3D STED) with two channels","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia --project=. scripts/ercontacts.jl -i <data> -o <output>","category":"page"},{"location":"tutorial/#Processing-a-dataset-with-multiple-channels","page":"Tutorial","title":"Processing a dataset with multiple channels","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia --project=. scripts/batch.jl -i <data> -o <output> -r \"*[1,2,3].tif\"","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This expects the data to be organized like so","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"- top folder\n  - replicate number\n    - cell type\n      - Seriesxyz \n        - ...0.tif\n        - ...1.tif\n        - ...2.tif \n        - ...3.tif\n        - ...4.tif","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In this case the \"*[1,2,3].tif\" parameters indicate that you want contacts between the pairs of channels of files ending with 1,2, and 3.  The code will check your data, and if those files are there, you will get output like","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"1--2/\n    replicate \n        cell \n            ...\n1--3/\n2--3/","category":"page"},{"location":"postprocessing/#Postprocessing","page":"Postprocessing","title":"Postprocessing","text":"","category":"section"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"Once the contact maps have been computed, you often need quantification and additional filtering.  For example, coverage, features descriptors, and so forth.","category":"page"},{"location":"postprocessing/#Aggregating-CSV-files","page":"Postprocessing","title":"Aggregating CSV files","text":"","category":"section"},{"location":"postprocessing/#Sampling-contacts","page":"Postprocessing","title":"Sampling contacts","text":"","category":"section"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"In scripts/runcubesamplingondataset.jl you'll find a script that samples contacts with a sliding window, to avoid long tail statistics dominating the conclusion of any analysis. The paper goes into more depth why this is beneficial.","category":"page"},{"location":"postprocessing/#Preprocessing-and-filtering","page":"Postprocessing","title":"Preprocessing and filtering","text":"","category":"section"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"The background filter removes ghost effects (bleedthrough). If you want to tune this without invoking the full pipeline, you can do so:","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"Suppose we want to filter all tif files ending with \"1.tif\" or \"2.tif\" , for z=1 to 1.1 in 0.25 steps, and then compute the object properties.","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"julia --project=.  scripts/segment.jl --inpath mydir -z 1.0 -Z 1.1 -s 0.25 -r \"*[1,2].tif\"","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"For each file, for each filtered value, it will generate:","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"mask.tif\nmasked.tif","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"For all the files, it will generate a CSV with columns, where each row is an object:","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"size (nr of non zero voxels)\nweighted (intensity sum of objects)\nminimum, Q1, mean, mediam, Q3, maximum, std, kurtosis : describes the intensity distribution of the object\nxyspan : major axis of 2D projection, pixels\nzrange : extent in Z\nzmidpoint : midpoint in Z slice\ndistancetocentroid: distance of this object's centroid to centroid of all objects, describes clustering\ndistancetocentroid_normalized: rescale the above to 0-1, where 0 is centroid of all objects, 1 is maximum dispersion\ncentroidchannel{x,y,z} : in pixel coordinates, records the centroid of all objects, this is the reference for the distance computation\ncentroidobject{x,y,z} : in pixel coordinates, records the object centroid\nfilename : the source tif file name\nz : the z value used\neig1-3: PCA eigenvalues, the can be used for shape descriptors\neig1-3normalized: eigenvalues rescaled to 1","category":"page"},{"location":"clustercomputing/#Cluster-Computing","page":"Cluster Usage","title":"Cluster Computing","text":"","category":"section"},{"location":"clustercomputing/","page":"Cluster Usage","title":"Cluster Usage","text":"This page documents the usage of the contact detection on [SLURM] based clusters.","category":"page"},{"location":"clustercomputing/","page":"Cluster Usage","title":"Cluster Usage","text":"Examples using the cluster at the Digital Alliance are provided here","category":"page"},{"location":"clustercomputing/","page":"Cluster Usage","title":"Cluster Usage","text":"Cluster specific scripts can also be found here","category":"page"},{"location":"installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"This project is developed using Julia.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"For ease of use and to maximize reproducibility we also provide container images using Singularity.","category":"page"},{"location":"installation/#Source-code","page":"Installation","title":"Source code","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"The fastest and easiest way is to clone the repository using Git. Alternatively you can download zipped releases. The below assumes you have Julia 1.9 or higher installed.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"In an empty, new folder:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"git clone https://github.com/bcardoen/SubPrecisionContactDetection.jl.git\ncd SubPrecisionContactDetection.jl\njulia --project=. -e `using Pkg; Pkg.build(); Pkg.test()`","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"This downloads the source code in a subfolder, builds it with all dependencies and tests it.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Once you have this, you can use either a terminal or an IDE (e.g. Visual Studio Code) to work with the source code to process new datasets.","category":"page"},{"location":"installation/#Singularity/Apptainer","page":"Installation","title":"Singularity/Apptainer","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"We provide a preconfigured container image with all dependencies here.","category":"page"},{"location":"#SubPrecisionContactDetection.jl-Documentation","page":"SubPrecisionContactDetection.jl Documentation","title":"SubPrecisionContactDetection.jl Documentation","text":"","category":"section"},{"location":"","page":"SubPrecisionContactDetection.jl Documentation","title":"SubPrecisionContactDetection.jl Documentation","text":"Welcome to the documentation for this package. Please see the sidebar for relevant sections.","category":"page"},{"location":"","page":"SubPrecisionContactDetection.jl Documentation","title":"SubPrecisionContactDetection.jl Documentation","text":"Depth = 5","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"<a name=\"faq\"></a>","category":"page"},{"location":"faq/#Troubleshooting-and-FAQ","page":"Help and FAQ","title":"Troubleshooting & FAQ","text":"","category":"section"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"If you have any issues, please create an issue.","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"Make sure to include:","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"include OS, Julia version\ndescription of steps to reproduce\nbe concise yet complete","category":"page"},{"location":"faq/#Can-I-change-the-singularity-image-?","page":"Help and FAQ","title":"Can I change the singularity image ?","text":"","category":"section"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"Yes, if you clone the repository, and are using Linux, you need to do 2 things","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"edit singularity_recipes/recipe.def\nexecute buildimage # Needs sudo","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"./buildimage.sh","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"This will rebuild the image, first checking out the latest version of the code.","category":"page"},{"location":"faq/#System-requirements","page":"Help and FAQ","title":"System requirements","text":"","category":"section"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"Expected RAM usage for images of sizes 500x500x20 ~ 5GB RAM, 2000x2000x70: ~ 50GB RAM, and so on. By default, all of JULIANUMTHREADS cores will be used to run in parallel. > 8 is overkill, so set to 4-8 at most:","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"export JULIA_NUM_THREADS=4","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"On desktops this is unlikely to be an issue, but on a cluster node with > 64 cores you will probably get a slowdown if you exceed 8-12 cores.","category":"page"},{"location":"faq/#I-cloned-the-repo-but-I-get-conflicts-during-the-installation-?","page":"Help and FAQ","title":"I cloned the repo but I get conflicts during the installation ?","text":"","category":"section"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"First, make sure you install and clone in a clean environment:","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"mdkir mydir\ncd mydir\njulia\njulia> ]\n(@v1.x) pkg> activate .\n(@v1.x) pkg> update","category":"page"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"Do not use Julia < 1.7, there's no guarantee that deprecated APIs will still work, and performance and user friendliness of the e.g. the package manager alone make 1.7 the ideal baseline.","category":"page"},{"location":"faq/#Memory-usage","page":"Help and FAQ","title":"Memory usage","text":"","category":"section"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"Current memory usage is higher than it strictly needs to be because we generate a lot of intermediate steps. In principle we could reduce usage by x2 or more, but it would come at the cost of debugging/interpretability.","category":"page"},{"location":"faq/#Installation-gives-errors-on-MacOs","page":"Help and FAQ","title":"Installation gives errors on MacOs","text":"","category":"section"},{"location":"faq/","page":"Help and FAQ","title":"Help and FAQ","text":"MacOS + Conda has a bug where a certificate error triggers a cascade of errors. The errors can be ignored, including the failing tests, this is an optional part of the module. When the bug in conda is resolved, this issue should be resolved as well.","category":"page"},{"location":"output/#Output","page":"Generated output","title":"Output","text":"","category":"section"},{"location":"output/","page":"Generated output","title":"Generated output","text":"A brief overview of the generated output follows here:","category":"page"},{"location":"output/","page":"Generated output","title":"Generated output","text":"Files are saved with the prefix of the original: If you have files, say, ABC1.tif, and ABC2.tif, then output will be saved ABC_xyz.tif, and so forth.","category":"page"},{"location":"output/","page":"Generated output","title":"Generated output","text":"PREFIX_confidence_map\n.tif: Each voxel holds the confidence (significance, p-value [0-1]) of the corresponding contact.\nPREFIX__channel_\n[1,2].tif: The filtered mitochondria and ER channel. Use this to inspect if the z-filter removed too much or too little\nPREFIX__pre_split_raw\n.tif: Unfiltered contacts. At voxel x, y, z this will have a correlation value, and its significance can be found in confidence_map[x,y,z].\nPREFIX_pre_split_gradient\n.tif: Contacts with the gradient filter applied.\nPREFIX_split_eroded\n.tif: Erodes singular voxels that are below the precision of the system\nPREFIX_3_eroded_volumes_nonsplit\n.csv : Features computed on the contacts\nPREFIX_C1_objects\n.csv: A CSV file where each row describe the features of the segmented objects of that channel (1). So if the code ran on 01.tif and 02.tif, C1 will map to 01.tif, filtered, then processed.\nPREFIX_C2_objects\n.csv: A CSV file where each row describe the features of the segmented objects of that channel (1). So if the code ran on 01.tif and 02.tif, C1 will map to 01.tif, filtered, then processed.","category":"page"},{"location":"output/","page":"Generated output","title":"Generated output","text":"See also the postprocessing for further output. The remainder are debugging outputs that can be traced to their source code in the script.","category":"page"}]
}
