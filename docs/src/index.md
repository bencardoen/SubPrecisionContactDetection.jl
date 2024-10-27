# SubPrecisionContactDetection.jl Documentation

Welcome to the documentation for this package.
Please see the sidebar for relevant sections.


## Index

```@index
```

## Notes

### Building the documentation
```julia
julia -e 'using Pkg; Pkg.activate(); push!(LOAD_PATH, pwd());'
julia --project=docs/ -e 'using Pkg; Pkg.activate();  push!(LOAD_PATH,pwd());'

julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate();'
julia --project=docs/ --color=yes docs/make.jl
```