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
using ArgParse
import ImageMagick
import Logging
import Random
using ImageFiltering
using Distributions
using Base.Threads
import JLD2

using LoggingExtras, Dates


# Code reused from DataCurator.jl



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--inpath", "-i"
            help = "input folder"
            arg_type = String
            required = true
        "--inregex", "-r"
            help = "regex tiff files in inpath. Defaults to *[1,2].tif , so will pick up tif files endings with 1 or 2.tif"
            arg_type = String
            default = "*[1,2].tif"
        "--startz", "-z"
            help = "Start Z value, default = 1.0"
            arg_type = Float64
            default = 1.0
        "--endz", "-Z"
            help = "End Z value, default = 3.5"
            arg_type = Float64
            default = 3.5
        "--step", "-s"
            help = "Step, default = 0.1"
            arg_type = Float64
            default = 0.1
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
    
    inpath = parsed_args["inpath"]
    z = parsed_args["startz"]
    Z = parsed_args["endz"]
    s = parsed_args["step"]
    pattern = parsed_args["inregex"]
    if z <= 0
        @error "z should be > 0"
        return
    end
    if z > Z
        @error "z > Z !!"
        return
    end
    if s <= 0
        @error "Step of <= 0 is nonsense"
        return
    end
    filter_mcsdetect(inpath, z, s, Z, pattern)
    @info "Because the glass is already broken, it is more enjoyed -- Ajahn Chah"
end

runc()