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

function runc()
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
    # Check the glob of input and how many files it finds 
    inpath = parse_args["inpath"]
    files = Glob.glob(parse_args["channels"], inpath)
    if length(files) < 2 || length(files) > 8
        @error "Files is < 2 or > 8, please check your channel specification."
        return
    end
    ends = endings(files)
    cs, cis = combines(ends)
    for ((xi, yi), (x, y)) in zip(cis, cs)
        @info "Starting combination $x - $y"
        op = joinpath(parse_args["outpath"], "$(x)-$(y)")
        pa = copy(parsed_args)
        pa["outpath"] = op
        tiffiles = [files[xi], files[yi]]
        two_channel_contacts(pa, tiffiles)
    end
end

