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
try
    @info "Importing ..."
    PyCall.pyimport("kneed");
    PyCall.pyimport("skimage");
    @info "Done"
catch e
    @warn "Failed import $e -- installing"
    println("Failed import $e -- installing")
    Conda.pip_interop(true)
    Conda.pip("install", "kneed")
    Conda.pip("install", "scikit-image")
end
@info "Success!"
