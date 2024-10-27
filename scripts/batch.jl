# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# Copyright 2021-2024, Ben Cardoen
using Pkg; Pkg.activate(".")
using SubPrecisionContactDetection
using Images, Colors, DataFrames, CSV, Statistics, LinearAlgebra
import Glob
using ArgParse
import ImageMagick
import Logging
import Random
using ImageFiltering
using Distributions
using Base.Threads
import JLD2

using LoggingExtras, Dates

inpath="/home/bcardoen/cedar_data/test/MS_2024_09_19_TOM20_KDEL_PLIN_HUH7_OLEIC"
# readdir(inpath)
regex="*[0,1].tif"
prefixes, regexes = buildregex(inpath, regex)

t = mktempdir()
subd = joinpath(t, "t2")
mkpath(subd)
Images.save(joinpath(subd, "a01.tif"), rand(200, 200))
Images.save(joinpath(subd, "b02.tif"), rand(200, 200))
Images.save(joinpath(subd, "c03.tif"), rand(200, 200))
rx = "*[1,2].tif"
paths, regexes = buildregex(subd, rx)
paths == ["1--2"]
regexes == ["*[1,2].tif"]

# function buildregex(path::AbstractString, regex::AbstractString)
#     @debug "Building regex combinations for $(regex) in $(path)"
#     cs, _ = test_multichannel(path, regex)
#     ps = Vector{String}()
#     rs = Vector{String}()
#     @info "Have $(length(cs)) combos : $(cs)"
#     for csx in cs
#         out = "$(csx[1])--$(csx[2])"
#         r = "*[$(csx[1]),$(csx[2])].tif"
#         push!(rs, r)
#         push!(ps, out)
#     end
#     @debug "Created path postfixes $(ps)  and regexes $(rs)"
#     return ps, rs
# end

# function test_multichannel(path, regex)
#     fs = recursive_glob(regex, path)
#     prefix_dict = Dict()
#     for f in fs
#         d = dirname(f)
#         if d in keys(prefix_dict)
#             push!(prefix_dict[d], f)
#         else
#             prefix_dict[d] = [f]
#         end
#     end
#     ## Check that the prefix dict has the same files for each key
#     ks = keys(prefix_dict) |> collect
#     vn = [length(prefix_dict[k]) for k in ks]
#     if !check_equal(vn)
#         @error "Unequal number of matches!!!"
#         @error [prefix_dict[k] for k in ks]
#         throw(ArgumentError("Expecting same amount of files for each subdirectory, aborting."))
#     end
#     files = prefix_dict[ks[1]]
#     ends = endings(files)
#     @debug "File endings $(ends)"
#     cs, cis = combines(ends)
#     @debug "Found combos $(cs)"
#     return cs, cis
# end

# function check_equal(xs)
#     if length(xs) <= 2
#         return true
#     end
#     return all(xs[1] .== xs)
# end

# function run_script()
#     date_format = "yyyy-mm-dd HH:MM:SS"
#     timestamp_logger(logger) = TransformerLogger(logger) do log
#       merge(log, (; message = "$(Dates.format(now(), date_format)) $(basename(log.file)):$(log.line): $(log.message)"))
#     end
#     ConsoleLogger(stdout, Logging.Info) |> timestamp_logger |> global_logger
#     parsed_args = parse_commandline_ercontacts()
#     @info "Parsed arguments:"
#     for (arg,val) in parsed_args
#         @info "  $arg  =>  $val"
#     end
#     inpath = parsed_args["inpath"]
#     op = parsed_args["outpath"]
#     # Check if we're multichannel
#     # If so, create output with multichannel
#     # Recursive glob of the regex
#     # Sort by directory
#     # Get the first directory
#     # create the patterns
#     for replicate in readdir(inpath; join=true)
# 	    r = basename(replicate)
# 	    for celltype in readdir(replicate; join=true)
# 	        ct = basename(celltype)
# 	        for cell in readdir(celltype; join=true)
# 	            snr = basename(cell)
#                 pa = copy(parsed_args)
#                 # op = parsed_args["outpath"]
#                 opx = joinpath(op, r, ct, snr)
#                 if isdir(opx)
#                     @warn "WARNING: Output directory exists --> PLEASE CHECK if this intended."
#                 end
#                 mkpath(opx)
#                 @info "Creating output path $(opx)"
#                 pa["inpath"]=cell
#                 pa["outpath"]=opx
# 				two_channel_contacts(pa)
# 	        end
# 	    end
# 	end
# end

function test_multichannel(inpath, regex)
    fs = recursive_glob(inpath, regex)
    prefix_dict = Dict()
    for f in fs
        d = dirname(f)
        if d in keys(prefix_dict)
            push!(prefix_dict[d]
        else
            prefix_dict[d] = [f]
        end
    end
    k = keys(prefix_dict) |> collect
    files = k[1]
    if length(files) < 2
        @error "Unexpected nr of files, for multichannel you should have at a minimum 2, got $(files)"
    end
    ends = endings(files)
    cs, cis = combines(ends)
    return cs, cis
end

# run_script()
