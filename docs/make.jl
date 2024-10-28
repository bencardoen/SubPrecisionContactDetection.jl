using Documenter
push!(LOAD_PATH,"../src/")
makedocs(sitename="SubPrecisionContactDetection Documentation",  pages=[ "Tutorial" => "tutorial.md", "Parameter selection" => "parameters.md", "Cluster Usage" => "clustercomputing.md", "Installation" => "installation.md"])

deploydocs(repo = "github.com/bencardoen/SubPrecisionContactDetection.jl.git")