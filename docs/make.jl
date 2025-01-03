using Documenter
push!(LOAD_PATH,"../src/")
makedocs(sitename="SubPrecisionContactDetection Documentation",  pages=[ "Home" => "index.md", "Tutorial" => "tutorial.md", "Parameter selection and tuning" => "parameters.md", "Generated output" => "output.md", 
"Cluster Usage" => "clustercomputing.md", "Installation" => "installation.md", "Validation and Limitations" => "validation.md",
"Postprocessing" => "postprocessing.md", "Help and FAQ" => "faq.md"])

deploydocs(repo = "github.com/bencardoen/SubPrecisionContactDetection.jl.git")