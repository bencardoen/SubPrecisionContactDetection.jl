# Validation and Limitations

## Validation
How do you know the detected interaction is `real`?
Here we describe methods to validate your results.

### Ground truth
#### Voxel level
If you have voxel level ground truth, you're unlikely to be needing the contact detection anyway. 
This is the easiest, but most rare case.
#### Image level
Here you have weak supervision, for example you know that in cell A there is a difference in interaction.
You can then optimize the detection's [parameters](https://bencardoen.github.io/SubPrecisionContactDetection.jl/dev/parameters/) to maximize that difference.

### Simulation
If you know what biological interaction you wish to observe and under what modality, you can simulate this and test if the detection can detect the simulated interaction, and to what extent. 

### Stability
Stability of a function ``f`` is defined as 
``\vert f(x + \epsilon) - f(x) \vert \leq \delta`` for ``\epsilon, \delta \geq 0, s = \vert \epsilon -\delta \vert``. 

If ``s`` is small for all values, then you have a stable function. 
For example, stability can be with respect to parameters or noise.
Stability matters for validation because you do not want to publish results that only exist for ``\alpha=0.003`` and w=2, and no other values.
Ideally on representative (median) cells, you do a parameter sweep to show a consistent difference or a consistent recall/precision with respect to ground truth.

!!! note "Stability"
    We made the simplification of using addition here, that is rarely true, but the principle holds. In addition, you do not want the stability to be small, but predictable. For example, as ``\epsilon`` increases, how does ``\delta`` increase? If there is a limiting function (linear, quadratic, ...), then that too is consistent. If this random, that is the true worst case.

### Phantoms
You can test the detection on phantoms, e.g. physical or biological induced changes (e.g. SPLICS), where you alter the biology in such a way that forces organelles or proteins to interact.

## Limitations
There are several factors affecting the outcome of the algorithm, we briefly describe each of them with mitigations:
For a description of confounding factors in microscopy in the context for interaction, see this [paper](https://zenodo.org/records/14009143) where we describe the factors and list methods to resolve them.

### Signal to noise ratio
In microscopy noise is a complex, largely unknown, and non-additive perturbation to the image.
#### Effect
Noise, especially noise that introduces pixellation effects, will disrupt the detection.
This will be first visible in the ``\alpha`` filtered values, as SNR goes from 2 to 1, say, you may need to increase ``\alpha`` (significance).

#### Mitigation
Apply deconvolution, ideally in combination with empirical resolution estimation tools such as [FRC/FSC](https://imagej.net/plugins/fourier-ring-correlation).
You would pick the deconvolution parameters based on the smallest average FRC.

### Registration
Objects can move between scanlines or acquisitions, to correct this you need registration.

#### Effect
There is no way the contact detection can guess that registration is needed, this is modality specific. 

#### Mitigation
Apply registration, but be _very_ careful. Registration algorithms can assume the registered objects should perfectly align, this destroys the interaction.

### Resolution
The resolution limits what can be observed within one voxel, in one channel. 
#### Effect
Contact detection uses spatial-temporal context to improve on what can be detected, but there are limits.
If you have an empirical resolution of say 250nm, and no other confounding factors.

#### Mitigation
First, measure using e.g. FRC. Note that FRC tends to predict resolutions a bit too optimistic.
Second, use deconvolution or preprocessing to improve the SNR and resolution.
Third, validate your results.

!!! warn "Empirical != Theoretical"
    Empirical resulotion is the resolution defined by your actual data and acquisition, theoretical is a best, worst, or average case based on specific calibrated issues, assumptions, and theoretical models from the microscope.

### Statistics
If you acquire interactions of two organelles across a whole 3D cell, be forewarned that this will likely be very imbalanced data.
While the human visual system is _very_ good at focusing on entropy rich differences, statistics is frequency driven (mostly).
#### Mitigation
Carefully frame your research question to take into account imbalance. 
You can long-tail robust modelling, stratification, or extreme value theory models to account for this.
A simple example is looking only at the QXX quantiles of a measure to adaptively isolate the long tail of a distribution, for exmaple.