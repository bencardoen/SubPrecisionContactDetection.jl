using Pkg;
using Logging;
@info "Initiating build"
## We want the Conda local Python env, anything else is out of control
if !haskey(ENV, "PYTHON")
    @info "Python not set, using PyCall"
    ENV["PYTHON"] = ""
else
    @info "Python set to $(ENV["PYTHON"])"
end
Pkg.add("Conda")
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.build("Conda")
using PyCall
using Conda
# Conda.add("gcc=12.1.0"; channel="conda-forge")
# Conda.add("kneed"; channel="conda-forge")
# Conda.add("scikit-image")
# Conda.add("scipy=1.8.0")
# PyCall.pyimport("kneed");
# PyCall.pyimport("skimage");
try
    @info "Importing ..."
    PyCall.pyimport("kneed");
    PyCall.pyimport("skimage");
    @info "Done"
catch e
    @warn "Failed import $e -- installing"
    println("Failed import $e -- installing")
    Conda.pip_interop(true)
    Conda.add("gcc=12.1.0"; channel="conda-forge")
    #Pin this version, to avoid clashes with libgcc.34
    Conda.add("scipy=1.8.0")
    Conda.pip("install", "kneed")
    Conda.pip("install", "scikit-image")
end
@info "Success!"
