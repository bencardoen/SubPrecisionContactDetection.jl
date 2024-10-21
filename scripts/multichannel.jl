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



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--inpath", "-i"
            help = "input folder"
            arg_type = String
            required = true
        "--normalize"
            help = "If set, in --decon mode use channel normalization after background removal."
            action = :store_true
            default = false
            required = false
        "--cube-vesicle-size-ln"
            help = "The cube sampling analysis drops any contact adjacent to a mitochondria of size (ln) <= 9"
            arg_type = Int
            default = 9
            required = false
        "--cube-vesicle-sample-size"
            help = "The cube sampling analysis window size = 5*K x 5*K x K, with default K=5"
            arg_type = Int
            default = 5
            required = false
        "--cube-vesicle-intensity-mean"
            help = "The cube sampling analysis drops any contact adjacent to a mitochondria of mean intensity <= 0.2"
            arg_type = Float64
            default = 0.2
            required = false
        "--filtermode", "-f"
            help = "filtermode, one of ['arithmetic' (default), 'curvature', 'geometric', 'autotune']. For autotune -z is computed, so ignored if set. See --prc to scale autotuning"
            arg_type = String
            required = false
            default = "arithmetic"
        "--prc", "-p"
            help = "Precision recall balance. When filtermode == autotune, value > 1 will favor recall, < 1 precision, default 1"
            arg_type = Float64
            required = false
            default = 1.0
        "--denominator", "-e"
            help = "When filtermode == curvature, value > 1 increases recall, < 1 increases precision"
            arg_type = Float64
            required = false
            default = 1.0
        "--weighted"
            help = "Use intensity sum weighted filter for option volumethreshold, default false"
            action = :store_true
        "--sphericity"
            help = "Sphericity filter to detect vesicles"
            arg_type = Float64
            required = false
            default = -1.0
        "--nooutput", "-n"
            help = "If set, do not generate output (testing)"
            action = :store_true
        "--dry-run"
            help = "Only check arguments and directory access"
            action = :store_true
        "--save-numerical-data"
            help = "Save results not just in tiff/csv, but also in compressed binary .jld format"
            action = :store_true
        "--channels", "-c"
            help = "regex matching exactly 2 tiff files in inpath"
            arg_type = String
            default = "*[1,2].tif"
        "--outpath", "-o"
            help = "output folder"
            arg_type = String
            required = true
        "--mode"
            help = "mode, one of decon, non-decon, both, default: non-decon. If 'both', expects to have 4 tif files, 2 x channel, ch0x_R.tif, ch0x_D.tif"
            arg_type = String
            required = false
        "--deconvolved"
            help = "Set to true if input is deconvolved (default false). If !true, supply --sigmas 'x-y-z'"
            action = :store_true
        "--sigmas"
            help = "String of 3 floating point values 'x-y-z', e.g.1-2.5-3, of the σ of the PSF in x,y,z. Necessary if deconvolved=false, default 1-1-1"
            arg_type = String
            required = false
            default = "1-1-1"
        "--lpsigmas"
            help = "String of 3 floating point values 'x-y-z', e.g.1-2.5-3, to smooth the Laplacian operator if needed. Default 1-1-1. Do not change unless you know what you're doing."
            arg_type = String
            required = false
            default = "1-1-1"
        "--windowsize", "-w"
            help = "window size, stride of the response field. A (w*2+1)^dims window is used for the correlation. Larger strides are more expensive to compute."
            arg_type = Int
            default = 1
        "--radius"
            help = "Interpret volume threshold as radius of sphere, not volume (default false)"
            action = :store_true
        "--zscore", "-z"
            help = "z score to segment the 3D stack (μ + z σ), only applies if deconvolution is true"
            arg_type = Float64
            default = 3.0
        "--volumethreshold", "-v"
            help = "Drop segments where the volume (voxels) is < threshold"
            arg_type = Int
            default = 0
        "--volumethresholdchannel", "-c"
            help = "Channel to apply threshold to, 1(default), 2, or 3(both)"
            arg_type = Int
            default = 1
        "--dimension"
            help = "3(D) (default), or 2 (XY)"
            arg_type = Int
            default = 3
        "--minzslice", "-m"
            help = "Start Z dimension at this index"
            arg_type = Int
            default = 1
        "--alpha"
            help = "Minimum Type I error acceptable"
            arg_type = Float64
            default = 0.05
        "--beta"
            help = "Minimum Type II error acceptable. β = 1 - stat power"
            arg_type = Float64
            default = 0.05
        "--maxzslice", "-M"
            help = "Start Z dimension at this index"
            arg_type = Int
        "--skipstats", "-s"
            help = "If true, skip computing statistics / features on contacts"
            action = :store_true
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

