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
using ProgressMeter

using LoggingExtras, Dates

function run_script()
    date_format = "yyyy-mm-dd HH:MM:SS"
    timestamp_logger(logger) = TransformerLogger(logger) do log
      merge(log, (; message = "$(Dates.format(now(), date_format)) $(basename(log.file)):$(log.line): $(log.message)"))
    end
    ConsoleLogger(stdout, Logging.Info) |> timestamp_logger |> global_logger
    parsed_args = parse_commandline_ercontacts()
    @info "Parsed arguments:"
    for (arg,val) in parsed_args
        @info "  $arg  =>  $val"
    end
    inpath = parsed_args["inpath"]
    rootpath = parsed_args["outpath"]
    rootregex = parsed_args["inregex"]
    paths, regexes = buildregex(inpath, rootregex)
    @showprogress for (pat, reg) in zip(paths, regexes)
        @showprogress  for replicate in readdir(inpath; join=true)
            r = basename(replicate)
            @showprogress  for celltype in readdir(replicate; join=true)
                ct = basename(celltype)
                @showprogress  for cell in readdir(celltype; join=true)
                    # Configure output
                    snr = basename(cell)
                    pax = copy(parsed_args)
                    alpha = pax["alpha"]
                    opx = joinpath(rootpath, pat, r, ct, snr, "$(alpha)")
                    @info "Output will be saved in $(opx)"
                    if isdir(opx)
                        @warn "WARNING: Output directory exists --> PLEASE CHECK if this intended."
                    end
                    @info "Creating output path $(opx)"
                    mkpath(opx)
                    # Update the arguments
                    pax["inpath"]=cell
                    pax["outpath"]=opx
                    pax["inregex"] = reg
                    try
        			    two_channel_contacts(pax)
                    catch e
                        @error("Failed executing due to $(e)")
                    end
                end
            end
        end
    end
end

run_script()