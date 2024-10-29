using Documenter
push!(LOAD_PATH,"../src/")
makedocs(sitename="SubPrecisionContactDetection Documentation",  pages=[ "Tutorial" => "tutorial.md", "Parameter selection and tuning" => "parameters.md", "Generated output" => "output.md", 
"Cluster Usage" => "clustercomputing.md", "Installation" => "installation.md", "Postprocessing" => "Postprocessing.md", "Help and FAQ" => "faq.md"])

deploydocs(repo = "github.com/bencardoen/SubPrecisionContactDetection.jl.git")