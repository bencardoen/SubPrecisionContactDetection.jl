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
using SubPrecisionContactDetection
using Images, Colors, DataFrames, CSV, Statistics, LinearAlgebra
import Glob
using ProgressMeter
using ArgParse
import ImageMagick
import Logging
import Random
using ImageFiltering
using Distributions
using Base.Threads
import JLD2

using LoggingExtras, Dates


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--inpath", "-i"
            help = "input folder"
            arg_type = String
            required = true
        "--outpath", "-o"
            help = "output folder. If not used, output is saved in place. If used, output is collapsed into this path. DO NOT use if filenames can clash"
            arg_type = String
            default = ""
            required = false
        "--recursive"
            help = "Search in all subdirectories of input folder. CAUTION: Can generate extreme amounts of output. Default off."
            action = :store_true
            default = false
        "--inregex"
            help = "Pattern to match, default to *.tif "
            arg_type = String
            default = "*.tif"
        "--add-shape", "-s"
            help = "Eigenvalue based shape computation is expensive for 3D shapes. For large objects this should be disabled. Default false."
            action = :store_true
            default = false
    end

    return parse_args(s)
end

function runc()
    date_format = "yyyy-mm-dd HH:MM:SS"
    timestamp_logger(logger) = TransformerLogger(logger) do log
      merge(log, (; message = "$(Dates.format(now(), date_format)) $(basename(log.file)):$(log.line): $(log.message)"))
    end
    ConsoleLogger(stdout, Logging.Info) |> timestamp_logger |> global_logger
    parsed_args = parse_commandline()
    @info "Parsed arguments:"
    for (arg,val) in parsed_args
        @info "  $arg  =>  $val"
    end
    
    # inpath = parsed_args["inpath"]
    # pattern = parsed_args["inregex"]
    op = parsed_args["outpath"]
    if op != "" && (! isdir(op))
        @warn "Output path does not exist, creating ..."
        mkdir(op)
        @warn "... done"
    end
    fs = recursive_glob(parsed_args["inregex"], parsed_args["inpath"])
    @info "Found a total of $(length(fs)) files"
    @debug fs
    dfs = []
    @showprogress for f in fs
        i = Images.load(f)
        df = describe_objects(i)
        df[!,:filename] .= f
        push!(dfs, df)
    end
    DFX = vcat(dfs...)
    @info "Writing CSV with $(size(DFX)) rows and cols"
    # DFX = vcat(dfs...)
    CSV.write(joinpath(op,"stats.csv"), DFX)
    @info "... Done"
    @info "Because the glass is already broken, it is more enjoyed -- Ajahn Chah"
end

runc()