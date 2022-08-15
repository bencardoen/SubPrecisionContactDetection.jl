using Pkg;
using Logging;
@info "Initiating build"
## We want the Conda local Python env, anything else is out of control
ENV["PYTHON"] = ""
Pkg.add("Conda")
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.build("Conda")
using PyCall
using Conda
Conda.add("gcc=12.1.0"; channel="conda-forge")
Conda.add("kneed"; channel="conda-forge")
Conda.add("scikit-image")
Conda.add("scipy=1.8.0")
PyCall.pyimport("kneed");
PyCall.pyimport("skimage");
@info "Success!"
