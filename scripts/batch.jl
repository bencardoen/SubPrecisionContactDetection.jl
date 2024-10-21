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
    # TODO setup multichannel mode
    inpath = parsed_args["inpath"]
    op = parsed_args["outpath"]
    for replicate in readdir(inpath; join=true)
	    r = basename(replicate)
	    for celltype in readdir(replicate; join=true)
	        ct = basename(celltype)
	        for cell in readdir(celltype; join=true)
	            snr = basename(cell)
                pa = copy(parsed_args)
                # op = parsed_args["outpath"]
                opx = joinpath(op, r, ct, snr)
                if isdir(opx)
                    @warn "WARNING: Output directory exists --> PLEASE CHECK if this intended."
                end
                mkpath(opx)
                @info "Creating output path $(opx)"
                pa["inpath"]=cell
                pa["outpath"]=opx
				two_channel_contacts(pa)
	        end
	    end
	end
end


run_script()
