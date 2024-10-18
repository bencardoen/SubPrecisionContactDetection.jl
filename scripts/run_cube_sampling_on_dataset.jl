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
using ArgParse, SPECHT, SubPrecisionContactDetection, Images, CSV, Statistics, DataFrames, Glob,  ERGO

### CLI tool to run the sampled statistics on the contacts.

function selectct(cts, selector)
    for ct in cts
        if occursin(selector, ct)
            return ct
        end
    end
    @assert false
end



function extract_stats(tiffs, replicate, celltype, serienr, config)
    mito = tiffs[(celltype, serienr, replicate)]["mito"]
    contacts = tiffs[(celltype, serienr, replicate)]["contacts"]
    contact = selectct(contacts,"gradient")
    @info "Loading images"
    mimg, cimg = [Images.load(i) for i in [mito, contact]]
    @info "Filtering"
    out, ms, ls = filtermito(mimg, cimg, config["ls"], config["rmv"])
    zs = config["zs"]
    @info "Walking $(5*zs) x $(5*zs) x $zs over mito / ct"
    OM = ERGO.tomask(out)
    res = walk_cube(out, OM.*cimg, ERGO.aszero(cimg), zs*5, zs*5, zs, 1)
    df, me, ce = res
    nzdf = filter(:mtsurface => !iszero, df) # new data frame
    return me, ce, ms, ls, nzdf, out, mimg, cimg.*OM
end



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
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
        "--outpath", "-o"
            help = "output folder"
            arg_type = String
            required = true
		"--inpath", "-i"
            help = "input folder"
            arg_type = String
            required = true
    end

    return parse_args(s)
end


function run()
	@debug "Todo process both channels, not just mito."
	parsed_args = parse_commandline()
	inpath = parsed_args["inpath"]
	outpath = parsed_args["outpath"]
	cubesize= parsed_args["cube-vesicle-sample-size"]
	z=3
	alpha=0.05
	zc=cubesize
	vesiclesizeln= parsed_args["cube-vesicle-size-ln"]
	vesicleint = parsed_args["cube-vesicle-intensity-mean"]
	tiffs = Dict()
	for replicate in readdir(inpath; join=true)
	    r = basename(replicate)
	    for celltype in readdir(replicate; join=true)
	        ct = basename(celltype)
	        for cell in readdir(celltype; join=true)
	            snr = basename(cell)
				for alphaval in readdir(cell; join=true)
					_alpha = basename(alphaval)
					
					ai = tryparse(Float64, _alpha)
					if ai != alpha
						@warn "Skpping $ai"
					end
					c1 = Glob.glob("*channel_1.tif", alphaval)[1]
					c2 = Glob.glob("*channel_2.tif", alphaval)[1]
					cts = Glob.glob("*pre_split_*.tif", alphaval)
					@info "Replicate $r Cell $ct Serie $snr alpha=$ai"
					#@assert(false)
					tiffs[(ct, snr, r)] = Dict("mito" => c1, "er" => c2, "contacts" =>cts)
				end
	        end
	    end
	end
	@info "Have a total of $(length(keys(tiffs)|>collect)) cells"
	config = Dict("ls"=>vesiclesizeln, "rmv"=>vesicleint, "zs"=>cubesize)
	RD = []
	for (c, s, r) in keys(tiffs)
	    @info "Cell $c SNR $(s) Rep $(r)"
	    est = extract_stats(tiffs, r, c, s, config)
	    resultshts1 = est[5]
	    filtered_contacts, filtered_mito = est[8], est[6]
	    @info "Saving tifs"
	    Images.save(joinpath(outpath, "filtered_mito_celltype_$(c)_replicate_$(r)_serienr_$(s).tif"), filtered_mito)
	    Images.save(joinpath(outpath, "filtered_contacts_celltype_$(c)_replicate_$(r)_serienr_$(s).tif"), filtered_contacts)
	    resultshts1[!, "alpha"] .= alpha
	    resultshts1[!, "z"] .= z
	    resultshts1[!, "cube_z"] .= zc
	    resultshts1[!, "celltype"] .= c
	    resultshts1[!, "replicate"] .= r
	    resultshts1[!, "serienr"] .= s
	    push!(RD, resultshts1)
	end
	@info "Saving DF"
	DF = vcat(RD...)

	CSV.write(joinpath(outpath, "all.csv"), DF)
	@info "Done"
	## Save
end


run()
