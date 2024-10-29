# API

In order of complexity, we now list the key function calls that compose the pipeline.

## Full pipeline
The main function that handles the pipeline end-to-end:
```julia
two_channel_contacts(parsed_args, tiffiles=nothing)
```

See [Parameter selection and tuning.](@ref) for configuring the arguments.


## Filtering

### Gradient Filter
Due to pixellation you can induce shadowing effects. These have an exact mathematical precondition, so we can filter them out.

```julia
gradientfilter!(correlation, laplacian1 , laplacian2)
```
Let the 3rd derivative of the image X be ``\delta^3_X``. The filter is then defined as keeping the values where this equation is true:

```math
    \neq ((\vert \vert \delta^3_X \vert \vert  \times \vert \vert \delta^3_Y \vert \vert)==0 \land \vert \vert \delta^3_X \vert \vert  + \vert \vert \delta^3_Y \vert \vert)!=0
```

## Image correlation
```julia
sp2d(img1, img2, stride)
```

Returns a correlation array, z-values that can be processed into significance, and a stride (`windowsize`)

```julia
sp3d(img1, img2, stride)
```

These run multithreaded, using as many cores as Julia is allowed to use.