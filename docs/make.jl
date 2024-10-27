using Documenter, SubPrecisionContactDetection
push!(LOAD_PATH,"../src/")
makedocs(sitename="SubPrecisionContactDetection Documentation",  pages=["API" => "functions.md", 
"Tutorial" => "tutorial.md", "Parameter selection" => "parameters.md", "Cluster Usage" => "clustercomputing.md", "Installation" => "installation.md"])
