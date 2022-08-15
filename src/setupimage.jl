using SubPrecisionContactDetection
using Pkg
using PackageCompiler
Pkg.activate(".")
@info pwd()
create_sysimage(sysimage_path="sys_img.so", include_transitive_dependencies=false, cpu_target="generic", precompile_statements_file="dc_precompile.jl")
