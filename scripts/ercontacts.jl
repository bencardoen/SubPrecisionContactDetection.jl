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
# Copyright 2021-2022, Ben Cardoen
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
        "--inregex", "-r"
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


function run_script()
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
    # pass remainder to refactored function
    twochannelcontacts(parsed_args)
end

function twochannelcontacts(parsed_args, tiffiles=nothing)
    # date_format = "yyyy-mm-dd HH:MM:SS"
    # timestamp_logger(logger) = TransformerLogger(logger) do log
    #   merge(log, (; message = "$(Dates.format(now(), date_format)) $(basename(log.file)):$(log.line): $(log.message)"))
    # end
    # ConsoleLogger(stdout, Logging.Info) |> timestamp_logger |> global_logger
    # parsed_args = parse_commandline()
    # @info "Parsed arguments:"
    # for (arg,val) in parsed_args
        # @info "  $arg  =>  $val"
    # end
    dimension = parsed_args["dimension"]
    if dimension != 3
        if dimension != 2
            error("Invalid dimenions $dimension , should be 2, 3")
        end
        @warn "Using XY axis only, input expected to be 2-dimensional"
    end
    inpath = parsed_args["inpath"]
    mode = parsed_args["filtermode"]
    allowed = ["autotune", "geometric", "arithmetic", "curvature"]
    z = parsed_args["zscore"]
    prc = parsed_args["prc"]
    denom = parsed_args["denominator"]
    if !(mode in allowed)
        @error "$(mode) not a valid choice"
        exit(-1)
    else
        if mode != "autotune"
            if mode != "curvature"
                if z < 0.01
                @error "Invalid z score $(z)"
                exit(-1)
                end
                @info "Using $(mode) with z = $(z)"
            else
                if denom <= 0
                    @error "Invalid d value <= 0 $(z)"
                    exit(-1)
                end
                @info "Using $(mode) with denominator = $(denom)"
            end
        else
            if prc <= 0
                @error "Invalid prc $(prc)"
                exit(-1)
            end
            @info "Using autotune with PRC $(prc)"
        end
    end
    sigmas = [nothing, nothing, nothing]
    deconvolved = parsed_args["deconvolved"]
    if !deconvolved
        sigmas = parsesigmas(parsed_args["sigmas"])
    end
    lpsigmas = parsesigmas(parsed_args["lpsigmas"])
    outpath = parsed_args["outpath"]
    if ! isdir(outpath)
        @warn "$(outpath) does not exists, creating"
        mkpath(outpath)
    end
    w = parsed_args["windowsize"]
    volth = parsed_args["volumethreshold"]
    @info "Using volume threshold of $volth"
    if volth < 0
        @error "Negative volume threshold? $volth"
        return
    end
    if (w < 1)
        @error "Invalid window size $w should be >=1"
        exit(-1)
    end
    if (w > 3)
        @warn "Window of $w is equivalent to $((w*2+1)^3) -- expect cubic increase in runtime !!!"
    end
    vtch = parsed_args["volumethresholdchannel"]
    if  vtch > 3 || vtch < 0
        @error "Invalid parameter $(vtch)"
        return
    end
    stride = w
    @assert(stride >= 1)
    # Allow to be overridden by caller so multicontacts does the right thing
    if isnothing(tiffiles)
        tiffiles = Glob.glob(parsed_args["inregex"], inpath)
    else
        @info "Caller passed in exact files $tiffles"
    end
    @info "Found $(length(tiffiles)) in $inpath"
    @assert(length(tiffiles) == 2)
    @assert(tiffiles[1] < tiffiles[2])
    path_components = splitpath(tiffiles[1])
    cell_dir, filename, treatment_dir  = path_components[end-1], path_components[end], path_components[end-2]
    prefix = "$(treatment_dir)_$(cell_dir)_$(splitext(filename)[1])"
    @info "Using $(prefix) to save files..."
    if parsed_args["dry-run"] == true
        @info "Dry run complete, quitting"
        return
    end
    @info "Loading images"
    _im1, _im2 = Images.load(tiffiles[1]), Images.load(tiffiles[2])
    im1 = splitchannel(_im1)
    im2 = splitchannel(_im2)
    raw = nothing
    rawregex = "*[1,2]_raw.tif"
    if parsed_args["mode"] == "both"
        @assert false
        @error("Revise support for combined mode")
        @warn "Mix"
        rtiffiles = Glob.glob(rawregex, inpath)
        @assert(length(rtiffiles) == 2)
        @assert(rtiffiles[1] < rtiffiles[2])
        r1 = splitchannel(Images.load(rtiffiles[1]))
        r2 = splitchannel(Images.load(rtiffiles[2]))
        raw = [r1, r2]
    end
    normalize = parsed_args["normalize"]
    @info "Processing 1 = $(size(im1)) $(eltype(im1)) 2 = $(size(im2)) $(eltype(im2))"
    if dimension == 2
        rawcontacts, rawmkcontacts, filteredcontacts, gradientcontacts, img_1f, img_2f, sigmap = compute_contacts_2d(im1, im2, k=z, w=stride, deconvolved=deconvolved,
        sigmas=sigmas, geometric=(mode == "geometric"), autotune=(mode == "autotune"), prc=prc, curvature=(mode == "curvature"), denominator=denom,
        alpha=parsed_args["alpha"], beta=parsed_args["beta"], raw=raw, lpsigmas=lpsigmas, normalize=normalize)
    else
        rawcontacts, rawmkcontacts, filteredcontacts, gradientcontacts, img_1f, img_2f, sigmap = compute_contacts_3d(im1, im2, k=z, w=stride, deconvolved=deconvolved,
        sigmas=sigmas, geometric=(mode == "geometric"), autotune=(mode == "autotune"), prc=prc, curvature=(mode == "curvature"), denominator=denom,
        alpha=parsed_args["alpha"], beta=parsed_args["beta"], raw=raw, lpsigmas=lpsigmas, normalize=normalize)
    end
    if parsed_args["nooutput"] == true
        @info "Skipping output"
        return
    end
    @info "Saving images"
    mito = Images.N0f16.(img_1f)
    Images.save(joinpath(outpath,"$(prefix)_confidence_map.tif"), Images.N0f16.(sigmap))
    Images.save(joinpath(outpath,"$(prefix)_channel_1.tif"), mito)
    Images.save(joinpath(outpath,"$(prefix)_channel_2.tif"), Images.N0f16.(img_2f))
    Images.save(joinpath(outpath,"$(prefix)_pre_split_raw.tif"), Images.N0f16.(rawcontacts))
    Images.save(joinpath(outpath,"$(prefix)_pre_split_gradient.tif"), Images.N0f16.(gradientcontacts))
    @info "Saving features of objects in channels"
    df_c1 = describe_objects(Images.N0f16.(img_1f))
    CSV.write(joirpath(outpath, "$(prefix)_C1_objects.csv"), df_c1)
    df_c2 = describe_objects(Images.N0f16.(img_2f))
    CSV.write(joirpath(outpath, "$(prefix)_C2_objects.csv"), df_c2)

    rawcontacts, rawmkcontacts, filteredcontacts = nothing, nothing, nothing
    GC.gc()
    erodedcontacts = gerode(gradientcontacts)
    GQ = Images.N0f16.(gradientcontacts .* erodedcontacts)
    Images.save(joinpath(outpath,"$(prefix)_pre_split_eroded.tif"), GQ)

    ## If 2D stop here
    if dimension == 2
        @info "Dimension = 2 : Feature computation & filtering not implemented for dim != 3"
        _df = reportvolumes2D(GQ)
        CSV.write(joinpath(outpath, "$(prefix)_$(w)_2D_eroded_surfaces_nonsplit.csv"), _df)
        @info "Dimension = 2 : Advanced filtering/postprocessing not implemented for 2D."
        return
    end

    _df, skeleton = reportvolumes(GQ, sigmap; mito=mito)
    if isnothing(_df)
        @debug "No components to process, skipping ..."
    else
        CSV.write(joinpath(outpath, "$(prefix)_$(w)_3_eroded_volumes_nonsplit.csv"), _df)
        Images.save(joinpath(outpath, "$(prefix)_$(w)_skeleton_contacts.tif"), skeleton)
    end
    #
    if parsed_args["radius"]
            @info "Using radius $volth"
            volth = radius_to_volume(volth)
            @info "is volume $volth"
    end

    @info "Filtering with volume threshold $volth"
    @assert vtch == 1
    GC.gc()
    lower, higher, lowercontacts, highercontacts = filter_channels(img_1f, GQ, volth, weighted=parsed_args["weighted"], intensityfilter=false, sphericity_threshold=parsed_args["sphericity"])

    _ls = ["vesicle", "non-vesicle"]
    mt = [lower, higher]
    con = [lowercontacts, highercontacts]
    for (il, _l) in enumerate(_ls)
        contacts = Images.N0f16.(con[il] .* GQ)
        ch = Images.N0f16.(mt[il] .* img_1f)
        @info "Saving filtered channel $(il) / 2"
        Images.save(joinpath(outpath,"$(prefix)_filtered_channel_1_$(_l).tif"), ch)
        @info "Saving filtered contacts $(il) / 2"
        Images.save(joinpath(outpath,"$(prefix)_$(_ls[il])_contacts.tif"), contacts)
        @info "Saving statistics $(il) / 2"
        _df, _  = reportvolumes(contacts, sigmap; mito=mito)
        if isnothing(_df)
            @debug "No components to process, skipping ..."
        else
            CSV.write(joinpath(outpath, "$(prefix)_$(_ls[il])_contacts.csv"), _df)
        end
    end
    @info "Saving config"
    JLD2.jldsave(joinpath(outpath, "metadata.jld2"), true; metadata=parsed_args)
    @info "Because the glass is already broken, it is more enjoyed -- Ajahn Chah"
end

run_script()
