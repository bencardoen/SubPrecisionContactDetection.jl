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
# Copyright 2020-2022, Ben Cardoen
module SubPrecisionContactDetection
import Pkg
import PyCall
using Distributions
using LinearAlgebra
import HypothesisTests
import Images
using ArgParse
import ERGO
# import ImageView # Crashes headlesss mode
using Statistics
import CSV
import JSON
using Logging
import Combinatorics
import DataFrames
import Glob
import Random
using Base.Threads
using ImageFiltering
using Match
using JLD2
using StatsBase
using ProgressMeter
using Colocalization
import ImageMorphology
import SPECHT
using ImageContrastAdjustment

export toct, getbox, edge_stack, binarize, spcor, magnitudegradient3d, computecontacts, normalizemaxmin, buildregex,
computeintensitycorrelation, recursive_glob,
summarize_spots, findchannel, sp, get_defaults,
compute_edges, reportimagequality, dtd_to_field, c3, shape_component, filter_mcsdetect, 
process_contact_stack3d, sp2d,
makespheres, loadimages, sp3d, reportvolumes, reportvolumes2D, process_contact_stack, filter_k, offset, mcc, clampt, ratefilter, computesurfaces,
festimparam, retainnegative, snr, planar, smoothclamp, sphericity, computefeatures, anisotropy, compute_contact_slice, normimg,
savegld, imgpca, readgld, dtocent, quantify_adj_mito, filter_channels, getextent, masktoindices, z_to_significance, magnitudegradients, gradientmagnitude, gradientfilter!,
indicestomask!, scoremasks, compute_contacts_3d, compute_contacts_2d, splitchannel, cntimg, alphabeta!, compute_contacts_deconvolved_images_2d,
vpdf, iterativemedian, spearsubsample, fastgaussian2d, fastguasian3d, sphere3d, offsetindices, estimatecurvature, describe_objects,
spear, denoise, volume_to_radius, radius_to_volume, qnorm, compute_min_r_for_sample_corr, get_components_diag, get_components_ndiag, get_components,
computesphericity, filterintensity, gerode, to3RGB, getskeleton, dropleq, randomcolorarray, getrandomcolor, minnz, normcontacts, getci,
gradientmagnitude, indexofdispersion, expandstack, reducestack, compute_sample_size_for_min_corr, gradientmagnitude3d,
reducedim, expanddim, parsesigmas, edge_from_border,
normalize_channel, filtermito, walk_cube, clipr, normalize_linear, combines, endings, two_channel_contacts, parse_commandline_ercontacts


### Define for the entire module how components are computed --> diag means touching any pixel
get_components_diag = mask -> Images.label_components(mask, length(size(mask))==2 ? trues(3,3) : trues(3,3,3))
get_components_ndiag = mask -> Images.label_components(mask)
get_components = get_components_diag


function normalize_linear(img)
    T = eltype(img)
    I = Float64.(img)
	@assert length(SPECHT.nz(I)) > 0
    A, B = minimum(SPECHT.nz(I)), maximum(SPECHT.nz(I))
    @debug "Normalizing($A-> $B to 0,1)"
    _i = copy(I)
    _i[_i .> 0] = (_i[_i .> 0] .- A) / (B-A)
    _i[_i .> 1] .= 1
    return T.(_i)
end

"""
	walk_cube(c1, ct, c3, stridex, stridey, stridez)

	Walks across c1/ct in x*y*z non overlapping cubes, recording content in dataframe.

	returns dataframe, borderc1, borderct
"""
function walk_cube(c1, ct, c3, stridex, stridey, stridez, border=1)
    X,Y,Z = size(c1)
	@debug X, Y, Z
    @assert(size(c1) == size(ct))
    @assert(size(c1) == size(c3))
    @assert(stridez > 0)
	dx, rx = divrem(X, stridex)
	dy, ry = divrem(Y, stridey)
	dz, rz = divrem(Z, stridez)
	@debug dx, dy, dz
	@debug "$X x $Y x $Z --> Total of $(dx*dy*dz) boxes"
	@debug rx, ry, rz
	c1border = edge_from_border(c1, border)
	borderct = Colocalization.tomask(ct .* c1border)
	if any([rx, ry, rz] .> 0)
		@warn " Slicing !"
	end

	df= DataFrames.DataFrame(mitvol = Int64[], mitosum=Float64[], contactvol=Int64[], contactsum=Float64[], ncontacts=Int64[], mtsurface=Int64[], ctsurface=Int64[])
	for _x in 1:dx
		for _y in 1:dy
			for _z in 1:dz
				c1n = c1[stridex*(_x-1)+1:stridex*_x, stridey*(_y-1)+1:stridey*_y, stridez*(_z-1)+1:stridez*_z]
				ctn = ct[stridex*(_x-1)+1:stridex*_x, stridey*(_y-1)+1:stridey*_y, stridez*(_z-1)+1:stridez*_z]
				ctbor = borderct[stridex*(_x-1)+1:stridex*_x, stridey*(_y-1)+1:stridey*_y, stridez*(_z-1)+1:stridez*_z]
				mtbor = c1border[stridex*(_x-1)+1:stridex*_x, stridey*(_y-1)+1:stridey*_y, stridez*(_z-1)+1:stridez*_z]
				mv, sm = length(c1n[c1n.>0]), sum(c1n)
				cv, sc = length(ctn[ctn.>0]), sum(ctn)
				surf_ct = length(ctbor[ctbor .> 0])
				surf_mt = length(mtbor[mtbor .> 0])
				@debug surf_ct surf_mt
				cnt = maximum(SubPrecisionContactDetection.get_components(Colocalization.tomask(ctn)))
				@debug [mv, sm, cv, sc, cnt]
				push!(df,[mv, sm, cv, sc, cnt, surf_mt, surf_ct])
			end
        end
    end
	return df, c1border, borderct
end

function clipr(img, seed=0, v=0.5)
	Random.seed!(seed)
	z = copy(img)
	X, Y, Z = size(img)
	rd = rand(X, Y, Z)
	z[rd .> v] .= 0
	return z
end

function erode_iter(img, N=1)
	M = Colocalization.tomask(img)
	em = copy(M)
	for _ in 1:N
		em = ImageMorphology.erode(em)
	end
	return em, M
end


"""
	edge_from_border(img, N=1)

	Creates a mask from img, then computes a boundary by eroding N times and remasking
	Returns a binary mask with 1 at border
"""
function edge_from_border(cellimg, N=1)
	eroded, mask = erode_iter(cellimg, N)
	result = Colocalization.aszero(mask)
	result[(eroded .== 0) .& (mask .== 1)] .= 1
	return result
end

"""
	expandstack(img, w)

	Given a 2D image, create a 3D stack of 2w+1 Z where each slice = img.
"""
function expandstack(img, w)
	S = size(img)
	@assert length(S) == 2
	@assert w > 0
	X, Y = S
	Z = 2*w + 1
	nimg = zeros(eltype(img), X, Y, 2*w+1)
	for z in 1:Z
		nimg[:,:,z] .= img
	end
	return nimg
end

"""
	reducestack(img, w)

	Get the w+1'th Z slice of img
"""
function reducestack(img, w)
	S = size(img)
	@assert(length(S) == 3)
	@assert w > 0
	Zs = 2*w + 1
	return copy(img[:,:,w+1])
end

function randomcolorarray()
    R, G, B = [Random.shuffle(X)./255 for _ in 1:3]
    return Images.RGB.(R, G, B)
end

"""
	getrandomcolor(seed)
	If seed != -1, set RNG with seed.
	Returns a random color in [0,1]^3
"""
function getrandomcolor(seed=-1)
    if seed != -1
        Random.seed!(seed)
    end
    return Images.RGB(rand(3)...)
end


function to3RGB(R, G, B)
    return SPECHT.tcolors([Images.N0f16.(R), Images.N0f16.(G), Images.N0f16.(B)])
end

function normcontacts(cts)
    CT = copy(cts)
    m = minnz(CT)
    R = 1 - m
    CT[CT.>0] .= (CT[CT .> 0] .- m)/R
    return CT
end

function minnz(xs)
    minimum(xs[xs .> 0])
end

function splitchannel(img)
    # if img is RGB, find the 1 non zero channel
    # if not exactly one channel is non zero, abort
    # if img is not RGB, return img
    if eltype(img) <: Images.RGB
		@info "RGB Image, trying to recover single channel"
		@debug "RGB image, looking for single NZ channel"
        CS = Images.channelview(img)
        FD = size(CS, 1)
        @assert(FD == 3)
        cnz = 0
        nz = 0
        for i in 1:FD
            zr = iszero(CS[i,:,:,:])
            if zr != true
                cnz += 1
                nz = i
            end
        end
        if cnz == 1
            @debug "RGB Image, selecting non zero channel"
            return copy(CS[nz, :, :, :])
        else
            @error "Image is RGB, with not exactly 1 non zero channel"
            error(-1)
        end
    else
        @debug "Grayscale images"
        return img
    end
end

"""
	compute_min_r_for_sample_corr(N, α, β)
	What is the minimum observable correlation r, given args.
	TODO: compute inverse, not approximate (but good enough for now)
"""
function compute_min_r_for_sample_corr(N, α, β)
    @assert N > 0
    @assert α > 0
    @assert β > 0
    # rs = 0:0.000001:1 |> collect
    # SN = r -> compute_sample_size_for_min_corr(r, α, β)
    # Ns = SN.(rs)
    # return rs[Ns .<= N][1]
	Z = sum(qnorm.([α, β]))
	Y = Z / (√(N-3))
	r = (exp(2*Y) - 1) / (exp(2*Y) + 1)
	return r
end


"""
	compute_sample_size_for_min_corr(r, α, β)
	What is the minimal sample size (Int) needed to register correlation r' >= r with type I α and type II β
"""
function compute_sample_size_for_min_corr(r, α, β)
	# Source https://www2.ccrb.cuhk.edu.hk/stat/other/correlation.htm
    ## Fisher's transform
    actual_z = 0.5 * log((1+r) / (1-r))
    # actual_z = atanh(r)
    za, zb = qnorm.([α, β])
    N = ((za + zb) / actual_z)^2 + 3
    return N
end

"""
	compute_contacts_3d
	Compute contacts between 2 3D volumes.
"""
function compute_contacts_3d(im1, im2; k=3, w=1, minz=nothing, maxz=nothing, sampleincrement=0, geometric=false,
	autotune=false, prc=1.0, curvature=false, denominator=1.0, deconvolved=false, sigmas=nothing, alpha=0.05, beta=0.05, raw=nothing, lpsigmas=[1,1,1], normalize=false)
	@assert(size(im1) == size(im2))
	@assert(length(size(im1))==3)
	@assert(isinteger(w))
	@assert(w>0)
	### Compute given a,b and N ((2w+1)^d) the min r
	minz = isnothing(minz) ? 1 : minz
	maxz = isnothing(maxz) ? size(im1, 3) : maxz
	im1 = im1[:,:,minz:maxz]
	im2 = im2[:,:,minz:maxz]
	if deconvolved
		@assert(k>=0)
		@debug "Deconvolution applied"
		if lpsigmas != [1, 1, 1]
			@warn "Laplacian sigmas not used in deconvolved images !!"
		end
		return compute_contacts_deconvolved_images_3d(im1, im2, w=w, k=k, geometric=false, autotune=false, prc=prc,
		curvature=curvature, denominator=denominator, alpha=alpha, beta=beta)
	else
		@debug "Deconvolution not applied"
		@assert(! isnothing(sigmas))
		return compute_contacts_nondeconvolved_images_3d(im1, im2, sigmas, w=w, alpha=alpha, beta=beta, raw=raw, lpsigmas=lpsigmas)
	end
end

function compute_contacts_2d(im1, im2; k=3, w=1, minz=nothing, maxz=nothing, sampleincrement=0, geometric=false,
	autotune=false, prc=1.0, curvature=false, denominator=1.0, deconvolved=false, sigmas=nothing, alpha=0.05, beta=0.05, raw=nothing, lpsigmas=[1,1,1], normalize=false)
	@assert(size(im1) == size(im2))
	@assert(length(size(im1))==2)
	@assert(isinteger(w))
	@assert(w>0)
	### Compute given a,b and N ((2w+1)^d) the min r
	# minz = isnothing(minz) ? 1 : minz
	# maxz = isnothing(maxz) ? size(im1, 3) : maxz
	# im1 = im1[:,:,minz:maxz]
	# im2 = im2[:,:,minz:maxz]
	if deconvolved
		@assert(k>=0)
		@debug "Deconvolution applied"
		if lpsigmas != [1, 1, 1]
			@warn "Laplacian sigmas not used in deconvolved images !!"
		end
		return compute_contacts_deconvolved_images_2d(im1, im2, w=w, k=k, geometric=false, autotune=false, prc=prc,
		curvature=curvature, denominator=denominator, alpha=alpha, beta=beta)
	else
		@assert false
		@debug "Deconvolution not applied"
		# @assert(! isnothing(sigmas))
		# return compute_contacts_nondeconvolved_images_3d(im1, im2, sigmas, w=w, alpha=alpha, beta=beta, raw=raw, lpsigmas=lpsigmas)
	end
end

"""
	parsesigmas(string)
	Decode a '-' separated string encoding of >= 0 floating point values into an array of 3 values.
	**note** exit(-1) on any failure
"""
function parsesigmas(strsig)
    splits = split(strsig, "-")
    if length(splits) != 3
        @error "Invalid preprocess option $(strsig), use format x-y-z"
        exit(-1)
    end
    sigmas = map(x -> tryparse(Float64, x), splits)
    if any(isnothing.(sigmas))
        @error "Invalid preprocess option $(strsig), use gaussian-x-y-z"
        exit(-1)
    else
        sx, sy, sz = sigmas
		if any(sigmas .<= 0)
			@error "∃ σ <= 0 !!! -- aborting"
        	exit(-1)
		end
        @debug "σ x $(sx) y $(sy) z $(sz)"
        return sigmas
    end
    exit(-1)
    return nothing
end

"""
	compute_contacts_nondeconvolved_images_3d
	For 2 given 3D volumes, raw, infer contacts.
	- w : window size
	- sigmas : PSF σ in x, y, z
	Uses a median filter to get rid of noise, Gaussian smoothing to replace deconvolution
"""
function compute_contacts_nondeconvolved_images_3d(img1, img2, sigmas; w=1, alpha=0.05, beta=0.05, raw=nothing, medianiters=2, lpsigmas=[1,1,1])
	@info "Preprocessing 0/2 -- median filter"
	## A first sweep to remove pixels that could be aggregated by the smoothing
	if isnothing(raw)
		@info "Only non decon images"
		S1N = iterativemedian(img1, 1, 1);
		S2N = iterativemedian(img2, 1, 1);
	else
		@warn "Using raw to mask out decon"
		S1N = iterativemedian(raw[1], 1, 1);
		S2N = iterativemedian(raw[2], 1, 1);
	end
	@info "Preprocessing 1/2"
	@debug sigmas
	S1NG = SubPrecisionContactDetection.denoise(Colocalization.tomask(S1N).*img1, 1,  "gaussian", sigmas);
	S2NG = SubPrecisionContactDetection.denoise(Colocalization.tomask(S2N).*img2, 1,  "gaussian", sigmas);
	@info "Preprocessing 2/2"
	S1NG[S1NG.<0].=0;
	S2NG[S2NG.<0].=0;
	@info "Differential 0/2"
	ll1 = imfilter(S1NG, Kernel.Laplacian((true, true, true)));
	ll2 = imfilter(S2NG, Kernel.Laplacian((true, true, true)));
	SPECHT.rlap!(ll1)
	SPECHT.rlap!(ll2)
	@info "Differential 1/2"
	@info "Using LP σ = $(lpsigmas)"
	grl1 = SubPrecisionContactDetection.denoise(ll1, 1, "gaussian", lpsigmas);
	grl2 = SubPrecisionContactDetection.denoise(ll2, 1, "gaussian", lpsigmas);
	@info "Differential Correlation"
	rawcontacts, zs, _ = sp3d(grl1, grl2, w)
	@info "Intensity Correlation"
	isf, izs, _ = sp3d(S1NG, S2NG, w)
	@info "Statistical Power"
	sigmap, ips  = z_to_significance(zs)[1], z_to_significance(izs)[1];
	ll1, ll2, zs, izs = nothing, nothing, nothing, nothing;
	@info "Reclaiming memory"
	GC.gc()
	@info "Filtering"
	SPECHT.rlap!(rawcontacts)
	SPECHT.rlap!(isf)
	# rawcontacts, isf = SPECHT.rlap(rawcontacts), SPECHT.rlap(isf)
	MS1 = iterativemedian(Images.N0f8.(S1NG), 1, 1);
	MS2 = iterativemedian(Images.N0f8.(S2NG), 1, 1);
	MSK = Colocalization.tomask(MS1.*MS2);
	filteredcontacts = copy(rawcontacts)
	isfx, filteredcontacts = alphabeta!(isf, ips, alpha, beta, w, 3),  alphabeta!(filteredcontacts, sigmap, alpha, beta, w, 3)
	filteredcontacts = MSK.*Colocalization.tomask(isfx).*filteredcontacts;
	gradientcontacts = gradientfilter!(filteredcontacts, grl1, grl2)
	@info "Completed --> PostProcessing"
	return rawcontacts.*MSK, rawcontacts .* MSK .* Colocalization.tomask(isf), filteredcontacts, gradientcontacts, S1NG.*Colocalization.tomask(MS1), S2NG.*Colocalization.tomask(MS2), sigmap
end


"""
	compute_contacts_deconvolved_images_3d
	For 2 given 3D volumes, deconvolved, infer contacts.
	- w : window size
	- k : segmentation step (roughly), ~3
	other options control segmentation
"""
function compute_contacts_deconvolved_images_3d(img1, img2; w=1, k=3, geometric=false, autotune=false, prc=1.0, curvature=false, denominator=1.0, alpha=.05, beta=0.05, normalize=false)
	@info "Filtering"
	img_1f, _ = filter_k(img1, k, geometric, autotune, 1.0, curvature, denominator)
	img_2f, _ = filter_k(img2, k, geometric, autotune, 1.0, curvature, denominator)
	I1 = copy(img1)
	I2 = copy(img2)
	if normalize
		@info "Using Normalization"
		I1 = normalize_linear(I1)
		I2 = normalize_linear(I2)
	end
	@info "Laplacian ..."
	ll1 = imfilter(I1, Kernel.Laplacian((true, true, true)));
	ll2 = imfilter(I2, Kernel.Laplacian((true, true, true)));
	ll1[ll1 .> 0] .= 0
	ll2[ll2 .> 0] .= 0
	rawcontacts, zs, _ = sp3d(ll1, ll2, w)
	isf, izs, _ = sp3d(img1, img2, w)
	sigmap, ips  = z_to_significance(zs)[1], z_to_significance(izs)[1];
	zs, izs = nothing, nothing;
	GC.gc()
	@debug "Use in place modifiers"
	# Use S1NG to mask
	rawcontacts, isf = SPECHT.rlap(rawcontacts), SPECHT.rlap(isf)
	rawcontacts[(img_1f .* img_2f) .== 0] .= zero(eltype(rawcontacts))
	filteredcontacts = copy(rawcontacts)
	isf, filteredcontacts = alphabeta!(isf, ips, alpha, beta, w, 3),  alphabeta!(filteredcontacts, sigmap, alpha, beta, w, 3)
	gradientcontacts = gradientfilter!(filteredcontacts, ll1, ll2)
	return rawcontacts, rawcontacts, filteredcontacts, gradientcontacts, img_1f, img_2f, sigmap
end

function ratefilter(GT, filtered)
    @assert size(GT) == size(filtered)
    z = zero(eltype(GT))
    on = oneunit(eltype(GT))
    FP = zeros(size(GT))
    FN = zeros(size(GT))
    TP = zeros(size(GT))
    TN = zeros(size(GT))
    FP[(filtered .> z) .& (GT .== z)] .= on
    TP[(filtered .> z) .& (GT .> z)] .= on
    FN[(filtered .== z) .& (GT .> z)] .= on
    TN[(filtered .== z) .& (GT .== z)] .= on
    total = reduce(*, size(GT))
	pos = sum(GT)
	@assert pos >= z
	neg = total - pos
	@assert neg >= z
    fp = sum(FP)/pos
    tp = sum(TP)/pos
    tn = sum(TN)/neg
    fn = sum(FN)/neg
    return FP, FN, TP, TN, fp, fn, tp, tn
end


function filtermito(mito, contacts, logsizemax, meanmitoint)
	mitomask = Colocalization.tomask(mito)
	out = copy(mito)
	contact_comps = SubPrecisionContactDetection.get_components(Colocalization.tomask(contacts))
	contact_indices = SubPrecisionContactDetection.get_components(contact_comps)[2:end]
	mito_comps = SubPrecisionContactDetection.get_components(mitomask)
	@debug "Have $(maximum(contact_comps)) contacts and $(maximum(mito_comps)) mito"
	indices = Images.component_indices(mito_comps)[2:end]
	logsize = []
	meanint = []
	for (i, ci) in enumerate(indices)
		mb = Float64.(mito[ci])
		sz = length(ci)
		sm = sum(mb)/sz
		_ls = log(sz)
		@debug sz, sm, _ls
		push!(meanint, sm)
		push!(logsize, _ls)
		if (_ls < logsizemax) && (sm < meanmitoint)
				@info " filtering "
				out[ci] .= 0
		end
	end
	return out, meanint, logsize
end

function normalize_channel(reference, other)
	@assert size(reference) == size(other)
	i1, i2 = adjust_histogram([reference, other], MidwayEqualization(nbins = 256))
	return i1, i2
end

function indexofdispersion(m, s)
	return (s)^2/m
end

"""
	Given a vector of values, map all values < -th to 0, all values > th to 1,
	and anything in between rescale to 0-1 range. (/2th)
"""
function clampt(_vec, th)
    vec = copy(_vec)
    vec .+= th
    vec[vec .< 0] .= 0.0
    vec[vec .> 2*th] .= 2*th
    vec ./= (2*th)
    return vec
end

function smoothclamp(_vec, k=64)
	vec = copy(_vec)
	return (1 .+ tanh.(vec .* k))./ 2
end

"""
	degrade(img, frombits, tobit)
	Simulates sparse encoding (e.g. 4-bit of 8-bit image)
"""
function degrade(img, from, to)
	res = copy(img)
	return Images.Normed{UInt, from}.(Images.Normed{UInt, to}.(res))
end

"""
	Return the probability that the values in zs (z-scores) are generated by a Normal(μ=0, σ=1)
	For convenience, return pdf, cdf and truncated pdf < alpha/2 (two tailed)
	Note Inf --> pdf(0), the mask is computed to omit Inf
"""
function z_to_significance(zs, alpha=0.05)
    pn = Distributions.Normal()
    pd, cd = pdf.(pn, zs), cdf.(pn, zs)
	msk = Colocalization.aszero(pd)
	msk[(pd .< (alpha/2)) .& (zs .!= Inf)].=1
	msk[zs .== Inf] .= 0
    return pd, cd, msk
end


function getci(xs, alpha)
    return HypothesisTests.confint(HypothesisTests.OneSampleTTest(xs), alpha, tail=:both)
end

"""
	gradientfilter
	Remove artifacts due to pixellation where F''' is 0 for either channel in either dim.
	Stores output in first argument
"""
function gradientfilter!(sf, ll1 , ll2)
    DIMS = length(size(sf))
    @assert((DIMS == 2) | (DIMS == 3))
    if length(size(sf)) == 2
		@info "Gradient filter with dim == 2"
        ## Compute magnitude of 3rd derivative
        δ2 = imgradients(ll2, KernelFactors.ando3);
        δ1 = imgradients(ll1, KernelFactors.ando3);
        mg1 = gradientmagnitude(δ1)
        mg2 = gradientmagnitude(δ2)
        ## If the magnitude of L' is zero for 1 channel, then we're not picking up a slope but pixellation
        # Where only 1 is zero, set to zero.
        sf[((mg2 .* mg1) .== 0) .& ((mg2 .+ mg1) .!= 0)] .= 0
        return sf
    else
		@info "Gradient filter with dim == 3"
        δ1x, δ1y, δ1z = imgradients(ll1, KernelFactors.ando3);
		δ2x, δ2y, δ2z = imgradients(ll2, KernelFactors.ando3);
        mg1 = magnitudegradients(δ1x, δ1y, δ1z)
        mg2 = magnitudegradients(δ2x, δ2y, δ2z)
        sf[((mg2 .* mg1) .== 0) .& ((mg2 .+ mg1) .!= 0)] .= 0
        return sf
    end
end


function compute_contacts_deconvolved_images_2d(img1, img2; w=1, k=3, geometric=false, autotune=false, prc=1.0, curvature=false, denominator=1.0, alpha=.05, beta=0.05, normalize=false)
	@info "Filtering"
	img_1f, _ = SubPrecisionContactDetection.filter_k(img1, k, geometric, autotune, 1.0, curvature, denominator)
	img_2f, _ = SubPrecisionContactDetection.filter_k(img2, k, geometric, autotune, 1.0, curvature, denominator)
	I1 = copy(img1)
	I2 = copy(img2)
	@assert length(size(img1)) == length(size(img2)) == 2
	if normalize
		@info "Using Normalization"
		I1 = normalize_linear(I1)
		I2 = normalize_linear(I2)
	end
	@info "Laplacian ..."
	ll1 = imfilter(I1, Kernel.Laplacian((true, true)));
	ll2 = imfilter(I2, Kernel.Laplacian((true, true)));
	ll1[ll1 .> 0] .= 0
	ll2[ll2 .> 0] .= 0
	rawcontacts, zs, _ = sp2d(ll1, ll2, w)
	isf, izs, _ = sp2d(img1, img2, w)
	sigmap, ips  = z_to_significance(zs)[1], z_to_significance(izs)[1];
	zs, izs = nothing, nothing;
	GC.gc()
	rawcontacts, isf = SPECHT.rlap(rawcontacts), SPECHT.rlap(isf)
	rawcontacts[(img_1f .* img_2f) .== 0] .= zero(eltype(rawcontacts))
	filteredcontacts = copy(rawcontacts)
	isf, filteredcontacts = alphabeta!(isf, ips, alpha, beta, w, 2),  alphabeta!(filteredcontacts, sigmap, alpha, beta, w, 2)
	gradientcontacts = gradientfilter!(filteredcontacts, ll1, ll2)
	return rawcontacts, rawcontacts, filteredcontacts, gradientcontacts, img_1f, img_2f, sigmap
end

function sp2d(img1, img2, stride)
    X, Y= size(img1)
    spx = zeros(Float32,X,Y)
	zs = zeros(Float32,X,Y)
    @assert(size(img1)==size(img2))
    @assert(stride >= 1)
    @assert(X>stride)
    @assert(Y>stride)
	N = (X-stride*2)*(Y-stride*2)
	# p = ProgressMeter.Progress(N, 0.5)
    @threads for x in 1+stride:X-stride
        for y in 1+stride:Y-stride
            @inbounds view1 = @view img1[x-stride:x+stride, y-stride:y+stride] # view saves 25% memory overhead
            @inbounds view2 = @view img2[x-stride:x+stride, y-stride:y+stride]
            _s, _z, _t = spear(view1, view2)
            @inbounds spx[x,y] = _s
			@inbounds zs[x,y] = _z
        end
    end
    return spx, zs, nothing
end

"""
	iterativemedian(image, iterations, window)
	Repeatedly apply image.*medianfilter(image)
"""
function iterativemedian(img, iters, w)
    cimg = copy(img)
    md = ones(eltype(img), size(img))
    for i in 1:iters
        md = SubPrecisionContactDetection.denoise(Colocalization.tomask(md).*cimg, w, "median",nothing)
    end
    return md
end

"""
	Return the magnitude of an array of gradients (grads[1] = array of gradients in dim 1
"""
function gradientmagnitude(grads)
    return sqrt.(mapreduce(x -> x.^2, (x,y) -> Float32.(x.+y), grads))
end

"""
	Map/reduce is too slow/mem intensive, explicit version for 3d
"""
function gradientmagnitude3d(dx, dy, dz)
    return sqrt.(dx.^2 .+ dy.^2 .+ dz.^2)
end


function scoremasks(G1, G2, G3, M1, M2, M3)
	res = Dict()
	GS = [G1, G2, G3]
	TS = [M1, M2, M3]
	for (i, (g, t)) in enumerate(zip(GS, TS))
		d, j = SPECHT.dice_jaccard(g, t)
		# Count
		cg = SubPrecisionContactDetection.cntimg(g)
		tg = SubPrecisionContactDetection.cntimg(t)
		res["$i"] = [d, j, cg, tg]
	end
	for (i,j) in Combinatorics.combinations(1:3, 2)
		gij = SPECHT.intersectimg(GS[i], GS[j])
		tij = SPECHT.intersectimg(TS[i], TS[j])
		d, k = SPECHT.dice_jaccard(gij, tij)
		# count
		cg = SubPrecisionContactDetection.cntimg(gij)
		ct = SubPrecisionContactDetection.cntimg(tij)
		res["$i$j"] = [d, k, cg, ct]
	end
	G123 = SPECHT.intersectimg(SPECHT.intersectimg(GS[1], GS[2]), GS[3])
	T123 = SPECHT.intersectimg(SPECHT.intersectimg(TS[1], TS[2]), TS[3])
	d, j = SPECHT.dice_jaccard(G123, T123, )
	cg, ct = SubPrecisionContactDetection.cntimg(G123), SubPrecisionContactDetection.cntimg(T123)
	res["123"] = [d, j, cg, ct]
	for k in keys(res)
		d, j, cg, ct = res[k]
		@debug "$k -> Dice $d Jaccard $j GT# $cg Mask# $ct"
	end
	return res
end

function cntimg(img)
	return maximum(get_components(img))
end

"""
	magnitudegradients
	Return the magnitude of the X/Y/Z gradients
"""
function magnitudegradients(gx, gy, gz=nothing)
    if isnothing(gz)
        return .√(gx .^ 2 .+ gy .^ 2)
    else
        return .√(gx .^ 2 .+ gy .^ 2 .+ gz .^ 2)
    end
end


"""
	spear(left, right)
	Return the spearman correlation between arguments, and the significance (z-test, t-test)
"""
function spear(view1, view2)
    N = length(view1)
    srt1 = sortperm(view1[:])
    srt2 = sortperm(view2[:])
    r = cor(srt1, srt2)
	z, t = computezcorr(r, N)
    return r, z, t
end

function computezcorr(r, N)
	# Fisher transform
	fr = atanh(r)
    z = sqrt((N-3)/1.06) * fr
    t = r * (sqrt((N-2)/(1-r^2)))
	return z, t
end

qnorm = α -> Distributions.quantile.(Distributions.Normal(0, 1), 1-(α/2))

function spearsubsample(view1, view2; increment=0)
	if increment==0 # Elide the copy
		return spear(view1, view2)
	else # Benchmarking shows this branch is compiled away, makes 0 difference.
		V1 = Images.imresize(view1, size(view1).+increment)
		V2 = Images.imresize(view2, size(view1).+increment)
		return spear(V1, V2)
	end
end

function offsetindices(indices, offset)
	nids = []
	for i in 1:length(indices)
		push!(nids, (indices[i] + offset))
	end
	return nids
end


function snr(intens)
	xs = Float64.(intens)
	return mean(xs)/std(xs)
end

function festimparam(xs, ys, x)
	@assert(length(xs) == length(ys))
	xd = x .- xs
	if all(xd .<= 0)
		@warn "x <= min table"
		return ys[1]
	end
	if all(xd .>= 0)
		@warn "x >= max table"
		return ys[end]
	end
	NX = length(xs)
	ND = length(xs[xd .< 0])
	i = NX - ND
	@debug "$i -> $(i+1)"
	xr = xs[i+1] - xs[i]
	@assert(xr > 0)
	ratio = (x - xs[i])/xr
	@debug ratio
	yr = ys[i+1] - ys[i]
	return ys[i] + ratio*yr
end

"""
	sphere3d
Return the x,y,z coordinates within radius distance of center
"""
function sphere3d(center, radius, X, Y, Z, skipradius)
	return [[x, y, z] for x in 1:X, y in 1:Y, z in 1:Z if skipradius <= dst([x, y, z], center) <= radius]
end

function dst(a, b)
	return sqrt(sum(( a.- b).^2))
end


"""
	vpdf(pdfs, coords)

Return sum of pdf on coord
"""
function vpdf(vs, coords)
	return SPECHT.vpdf(vs, coords)
end

"""
	fastgaussian2d(centers, covs, XN, YN)
Return an XN by YN array where each [x, y] = sum(pdf_gaus(center, cov)[x,y])
|centers| == |covs|, for each pair a Gaussian pdf is generated
"""
function fastgaussian2d(centers, covs, XN, YN)
	return SPECHT.fastgaussian2d(centers, covs, XN, YN)
end

"""
	fastgaussian3d(centers, covs, XN, YN, ZN)
Return an XN by YN array where each [x,y,z] = sum(pdf_gaus(center, cov)[x,y,Z])
|centers| == |covs|, for each pair a Gaussian pdf is generated
"""
function fastgaussian3d(centers, covs, XN, YN, ZN)
	vs =  [MvNormal(center, cv) for (center, cv) in zip(centers, covs)]
	Q = (vpdf(vs, [x, y, z]) for x in 1:XN, y in 1:YN, z in 1:ZN)
	return Q |> collect
end

function indicestomask!(indices, array)
	_n = oneunit(eltype(array))
	for ind in indices
		array[ind...] = _n
	end
	return array
end

function reportimagequality(imgfilename)
	img = Images.load(imgfilename)
	flattened = Float64.(img[:])
	nz = flattened[flattened .> 0]
	m, s, M = Statistics.mean(nz), Statistics.std(nz), Statistics.maximum(nz)
	iod = indexofdispersion(m, s)
	snr = m/s
	@debug "mean $(m) std $(s) max $(M) SNR $(snr) IOD $(iod)"
	return m, s, M, iod, snr
end

function dtd_to_field(λ, V)
    sa = length(size(V))
    @match sa begin
        4 => return dtd_to_field2(λ, V)
        5 => return dtd_to_field3(λ, V)
        _ => @assert(false)
    end
end

function dtd_to_field2(λ, V)
    X, Y, dims = size(λ)
    field = zeros(Float64, X, Y, 2)
    mg = zeros(Float64, X, Y)
    for y in 1:Y, x in 1:X
        λMi = argmax(λ[x,y,:])
        mg[x,y] = λ[x,y, λMi]
        field[x,y, :] = V[x,y,:, λMi]
    end
    return field, mg
end

function dtd_to_field3(λ, V)
    X, Y, Z, dims = size(λ)
    field = zeros(Float64, X, Y, Z, 3)
    mg = zeros(Float64, X, Y, Z)
    for z in 1:Z, y in 1:Y, x in 1:X
        λMi = argmax(λ[x,y,z,:])
        mg[x,y,z] = λ[x,y,z, λMi]
        field[x,y,z, :] = V[x,y,z,:, λMi]
    end
    return field, mg
end


function quantify_adj_mito(mitochannel, contactchannel)
	@assert size(mitochannel) == size(contactchannel)
	CM = Colocalization.tomask(contactchannel)
	skel = getskeleton(CM)
	c_comps = get_components(CM)
	c_inds = Images.component_indices(c_comps)[2:end]
	m_comps = get_components(Colocalization.tomask(mitochannel))
	m_inds = Images.component_indices(m_comps)[2:end]
	m_ls = Images.component_lengths(m_comps)[2:end]
	N = maximum(c_comps)
	if N == 0
		@warn "No contacts to process"
		return nothing, nothing
	end
	res = zeros(N, 3)
	for ic in 1:N
		# Get the mito on which contact is resting
		mitonrs = m_comps[c_inds[ic]]
		mnr = maximum(mitonrs)
		if mnr == 0
			@warn "Unattached contact $ic"
			MV, MW = 0.0, 0.0
		else
		# Get that mito's vol
			MV = m_ls[mnr]
			# And weighted
			MW = sum(mitochannel[m_inds[mnr]])
			@assert 0 < MW <= MV
		# Lookup skeleton area
		end
		SK = sum(skel[c_inds[ic]])
		@debug "Component $(ic) has skeleton area $SK on mito with volume $MV and weighted $MW"
		res[ic, :] .= MV, MW, SK
	end
	return res, skel
end

function denoise(img, window=1, type="median", sigmas=nothing)
    @match type begin
        "harmonic" => return denoisehm(img, window)
		"geometric" => return denoisegm(img, window)
		"median" => return denoisemd(img, window)
		"gaussian" => return denoisegs(img, sigmas)
        _ => @assert(false)
    end
end


function denoisegm(im, w)
	@debug "Geometric mean filter"
    mg= mapwindow(ERGO.gm, im, [w*2+1 for _ in 1:length(size(im))])
    mg[isnan.(mg)] .= zero(eltype(im))
	return mg
end


function denoisehm(im, w)
	@debug "Harmonic mean filter"
    mg= mapwindow(SPECHT.harmonicmean, im, [w*2+1 for _ in 1:length(size(im))])
    mg[isnan.(mg)] .= zero(eltype(im))
	return mg
end

function denoisemd(im, w)
	@debug "Median filter"
    mg= mapwindow(Statistics.median, im, [w*2+1 for _ in 1:length(size(im))])
    mg[isnan.(mg)] .= zero(eltype(im))
	return mg
end

function denoisegs(im, sigmas)
	@debug "Gaussian filter"
	@assert(length(sigmas) == length(size(im)))
	return imfilter(im, Kernel.gaussian((sigmas)))
end

function dropleq(mk, F)
    M = Colocalization.tomask(mk)
    CCS = get_components(M)
    l = Images.component_lengths(CCS)[2:end]
    li = Images.component_indices(CCS)[2:end]
    CF = copy(mk)
    for (i, ind) in enumerate(li)
        if l[i] < F
            M[ind] .= 0
        end
    end
    return M
end

function gerode(x)
	if length(size(x)) == 3
		@debug "Erosion with dim = 3"
	    res = copy(Colocalization.tomask(x))
	    a = falses(3,3,3)
	    a[1,2,2] = true
	    a[2,2,:] .= true
	    a[2,:,2] .= true
	    a[3,2,2] = true
	    return ImageMorphology.erode(res; dims=a)
	else
		@assert length(size(x))==2
		@debug "Erosion with dim == 2"
		return ImageMorphology.erode(Colocalization.tomask(x))
	end
end

"""
	process_contact_stack3d
	tiffiles : Array of strings to tiff images for channel 1, 2
	k : μ + k σ filter to correct for background bleeding
	w : filter width, 1 = 3x3x3, 2 = 5x5x5 etc
	minz, maxz = range of z slices
	sampleincrement : if > 0, how much to resample/interpolate to increase correlation filter size
"""
function process_contact_stack3d(tiffiles, k, w, minz=nothing, maxz=nothing; sampleincrement=0, geometric=false, autotune=false, prc=1.0, curvature=false, denominator=1.0, denoiseoption="off", densigmas=nothing)
    @assert(k > 0)
    @assert(w > 0)
	@assert(length(tiffiles) == 2)
	@error "DEPRECATED, compute_contacts_3d"
	@info "Processing with μ + $k σ, w $w , interpolation to $sampleincrement"
    @info "Loading ..."
    img_1 = Images.load(tiffiles[1])
    img_2 = Images.load(tiffiles[2])
	@info "Image dimensions = $(size(img_1))"
    @assert(size(img_1) == size(img_2))
    @assert(length(size(img_1)) == 3)
    _minz, _maxz = 1, size(img_1, 3)
    if !isnothing(minz)
        _minz = minz
    end
    if !isnothing(maxz)
        _maxz = maxz
    end
    @assert(1 <= _minz < _maxz <= size(img_1)[3])
    img_2 = img_2[:,:,_minz:_maxz]
    img_1 = img_1[:,:,_minz:_maxz]
    @assert(size(img_1) == size(img_2))
	if denoiseoption == "median"
		@debug "Denoising with option $denoiseoption $densigmas"
		# Median
		MS1 = iterativemedian(img_1, 1, 1)
		MS2 = iterativemedian(img_2, 1, 1)
		# Smooth
		img_1 = SubPrecisionContactDetection.denoise(Colocalization.tomask(MS1).*img_1, 1,  "gaussian", densigmas);
		img_2 = SubPrecisionContactDetection.denoise(Colocalization.tomask(MS2).*img_2, 1,  "gaussian", densigmas);
		@info "Laplacian ..."
		ll1 = imfilter(img_1, Kernel.Laplacian((true, true, true)));
    	ll2 = imfilter(img_2, Kernel.Laplacian((true, true, true)));
		rl1 = SPECHT.rlap(ll1)
		rl2 = SPECHT.rlap(ll2)
		ll1 = SubPrecisionContactDetection.denoise(rl1, 1, "gaussian", [1,1,1])
		ll2 = SubPrecisionContactDetection.denoise(rl2, 1, "gaussian", [1,1,1])
		img_1f = img_1
		img_2f = img_2
	else
		@info "Filtering ..."
		img_1f, th1 = filter_k(img_1, k, geometric, autotune, 1.0, curvature, denominator)
    	img_2f, th2 = filter_k(img_2, k, geometric, autotune, 1.0, curvature, denominator)
    	@assert(th1 > 0)
    	@assert(th2 > 0)
		@info "Laplacian ..."
		ll1 = imfilter(img_1, Kernel.Laplacian((true, true, true)));
    	ll2 = imfilter(img_2, Kernel.Laplacian((true, true, true)));
    	ll1[ll1 .> 0] .= 0
    	ll2[ll2 .> 0] .= 0
	end
    @info "Correlation ..."
    sf, zs, _ = sp3d(ll1, ll2, w, subsampleincrement=sampleincrement)
	@info "Intensity Correlation"
	fi, zis, _ = sp3d(img_1, img_2, w, subsampleincrement=sampleincrement)
	IC = SPECHT.rlap(fi)
	@info "Post-processing ..."
	rawsf = retainnegative(sf)
	sf = retainnegative(sf)
	if denoiseoption != "off"
		MS1 = iterativemedian(Images.N0f8.(img_1f), 1, 1)
		MS2 = iterativemedian(Images.N0f8.(img_2f), 1, 1)
		# Remove artifacts, intensity should be correlated, and neg cor
		qf = Colocalization.tomask(MS1.*MS2) .* Colocalization.tomask(IC) .* sf
	else
		sf[(img_1f .* img_2f) .== 0] .= zero(eltype(sf))
		qf = copy(sf)
		qf = gradientfilter!(qf, ll1, ll2)
	end
    return sf, img_1f, img_2f, zs, nothing, ll1, ll2, qf, rawsf
end

"""
	alphabeta!(values, p-values, α, β)
	values[values .< beta ad p-values > alpha] .= 0
	Used to filter an array of correlation values based on Type I/II error rates.
	TODO : fix b to its value
"""
function alphabeta!(array, ps, α, β, window=1, dim=3)
	## Compute R acceptable for ab at window
	@assert window >= 1
	@assert dim > 0
	@assert size(array) == size(ps)
	N = ((2*window)+1)^dim
	minr = compute_min_r_for_sample_corr(N, α, β)
	@debug "Minr $minr for $(α) $(β) with N $N"
	array[ps .>= α] .= 0
	array[array .< minr] .= 0
	return array
end

function retainnegative(array)
	ac = copy(array)
	z = zero(eltype(array))
	ac[ac .> z].= z
	ac[ac .< z].= z .- ac[ac .< z]
	return ac
end


function savegld(fname, rawcontacts, filteredcontacts, gradientcontacts, sigmap; tp=Float64)
	JLD2.jldsave(fname, true; rawcontacts=tp.(rawcontacts), filteredcontacts=tp.(filteredcontacts), gradientcontacts=tp.(gradientcontacts), sigmap=tp.(sigmap))
end


"""
	filter_mcsdetect : Filter a directory of tiff files with a sweep of z values
	channels is a regex to select which pattern of files to use.
	expected files should be 3D tiff 
	To use a single Z value, set start=stop.
	Object statistics are saved in a CSV file with the same name as the tiff file.
"""
function filter_mcsdetect(dir, start=1, step=0.1, stop=3, channels="*[0-2].tif", recursive=false, outpath="")
    @debug "Dir $dir sweep from $start → $stop in steps $step matching channels $channels in $recursive mode with $outpath"
    if recursive
		fs = recursive_glob(channels, dir)
	else
		fs = Glob.glob(channels, dir)
	end
    # @info fs
    @debug "Found $fs"
    for f in fs
        dfs = []
        i = Images.load(f)
        df = describe_objects(i)
        df[!,:z] .= NaN
        df[!,:filename] .= f
        push!(dfs, df)
        pt = splitpath(f)
        fn = pt[end]
        fne = splitext(fn)[1]
		savepath = pt[1:end-1]
		if outpath != ""
			savepath = splitpath(outpath)
		end
        for _z in start:step:stop
            fi, th = filter_k(i, _z)
            @debug "Threshold used for $(_z) : $(th)"
            m = bm(fi)
            @debug "Saving as mask_$(fne).tif) in $(joinpath(pt[1:end-1]...))"
            Images.save(joinpath(savepath...,"mask_$(_z)_$(fne).tif"), m)
            Images.save(joinpath(savepath...,"masked_$(_z)_$(fne).tif"), fi)
            df = describe_objects(fi)
            df[!,:z] .= _z
            df[!,:filename] .= f
            push!(dfs, df)
        end
        CSV.write(joinpath(savepath...,"stats_$(start)_$(step)_$(stop)_$(fne).csv"), vcat(dfs...))
    end
end


"""
	two_channel_contacts(args::Dict(String=>value), tiffiles=nothing)
	The main 2-channel contact function that combines preprocessing and postprocessing.
	See get_defaults for default arguments.
"""
function two_channel_contacts(parsed_args, tiffiles=nothing)
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
    if dimension == 3
        df_c1 = describe_objects(Images.N0f16.(img_1f))
        CSV.write(joinpath(outpath, "$(prefix)_C1_objects.csv"), df_c1)
        df_c2 = describe_objects(Images.N0f16.(img_2f))
        CSV.write(joinpath(outpath, "$(prefix)_C2_objects.csv"), df_c2)
    end 
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
    @info "Saving config"
    JLD2.jldsave(joinpath(outpath, "metadata.jld2"), true; metadata=parsed_args)
    @info "Because the glass is already broken, it is more enjoyed -- Ajahn Chah"
end


"""
	combines(xs)
	Returns a combination of the elements of xs and their indices
    Exmaple: 
    ```julia
    combs, inds = combines(["a", "b", "c"])
    ````
"""
function combines(xs)
    N = length(xs)
    inds = Combinatorics.combinations(1:N, 2)|> collect
    rs = []
    for ind in inds
        push!(rs, xs[ind])
    end
    return rs, inds
end


function parse_commandline_ercontacts()
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

function buildregex(path::AbstractString, regex::AbstractString)
    @debug "Building regex combinations for $(regex) in $(path)"
    cs, _ = test_multichannel(path, regex)
    ps = Vector{String}()
    rs = Vector{String}()
    @info "Have $(length(cs)) combos : $(cs)"
    for csx in cs
        out = "$(csx[1])--$(csx[2])"
        r = "*[$(csx[1]),$(csx[2])].tif"
        push!(rs, r)
        push!(ps, out)
    end
    @debug "Created path postfixes $(ps)  and regexes $(rs)"
    return ps, rs
end

function test_multichannel(path, regex)
    fs = recursive_glob(regex, path)
    prefix_dict = Dict()
    for f in fs
        d = dirname(f)
        if d in keys(prefix_dict)
            push!(prefix_dict[d], f)
        else
            prefix_dict[d] = [f]
        end
    end
    ## Check that the prefix dict has the same files for each key
    ks = keys(prefix_dict) |> collect
    vn = [length(prefix_dict[k]) for k in ks]
    if !check_equal(vn)
        @error "Unequal number of matches!!!"
        @error [prefix_dict[k] for k in ks]
        throw(ArgumentError("Expecting same amount of files for each subdirectory, aborting."))
    end
    files = prefix_dict[ks[1]]
    ends = endings(files)
    @debug "File endings $(ends)"
    cs, cis = combines(ends)
    @debug "Found combos $(cs)"
    return cs, cis
end

function check_equal(xs)
    if length(xs) <= 2
        return true
    end
    return all(xs[1] .== xs)
end



"""
	endings(fs)
	For a vector of filenames, return a vector of filename postfixes (integer). E.g. ["01.tif"] --> [1]
"""
function endings(fs)
    endings = []
    for f in fs
        f_split = splitpath(f)[end]
        f_name = splitext(f_split)[1]
        @debug "File name $(f_name)"
        if length(f_name) < 2
            @error "Unlikely filename $(f_name)"
            throw(ArgumentError("Unlikely filename $(f_name)"))
        end
        @debug "Detecting integer ending ...."
        intend = match(r"[0-9]+$", f_name)
        if isnothing(intend)
            throw(ArgumentError("Filename does not end with integer, can't make combos"))
        end
        ind = tryparse(Int, intend.match)
        if !isnothing(ind)
            push!(endings, ind)
        else
            throw(ArgumentError("Invalid channel name $(f_name)"))
        end
    end
    return endings
end


"""
    bm(xs)
    Return a binary (copy) mask of the array
"""
function bm(xs)
    ys = copy(xs)
    ys[ys .> 0] .= 1
    return ys
end

"""
    describe_objects(img::AbstractArray{T, 3}, shape=false) where {T<:Any}

    For a 3D array, describe objects and basic features
"""
function describe_objects(img::AbstractArray{T, 3}, shape=false) where {T<:Any}
    b = copy(img)
    b[b .> 0] .= 1
	## Changed 3-2 connectivity
	get_components_diag = mask -> Images.label_components(mask, length(size(mask))==2 ? trues(3,3) : trues(3,3,3))
    coms = get_components_diag(b)
    # lengths = Images.component_lengths(coms)[2:end]
    indices = Images.component_indices(coms)[2:end]
    boxes = Images.component_boxes(coms)[2:end]
    N = maximum(coms)
    w=zeros(N, 19)
	@info "Processing $N components"
	if N == 0
		@warn "NO COMPONENTS TO PROCESS"
		return nothing
	end
    @showprogress for ic in 1:N
        vals = img[indices[ic]]
		n = length(vals)
         # m, Q1, mx, med, Q3, M, std(ys), kurt = dimg(vals)
		w[ic, 2] = sum(vals)
		w[ic, 1] = n
		w[ic,3:10] .= _dimg(vals)
		w[ic, 11:13] .= getextent(boxes[ic])
        if shape
		    l1, l2, l3 = shape_component(coms, img, ic)	
		    w[ic, 14:16] .= l1, l2, l3
            if l1 > 0
                w[ic, 17:19] .= l1/l1, l2/l1, l3/l1
                @debug w[ic, 17:19]
            end
        end
		# if l1 == 0, l2, l3 are zero, the array is zero init, so they're already zero
	end
	columns = [:size, :weighted, :minimum, :Q1, :mean, :median, :Q3, :maximum, :std, :kurtosis, :xyspan, :zspan, :zmidpoint, :eig1, :eig2, :eig3, :eig1norm, :eig2norm, :eig3norm]
    df = DataFrames.DataFrame()
    for (i,c) in enumerate(columns)
        df[!,c] = w[:,i]
    end
    distances, nd, ctr, ctrs = _dtocent(coms)
    df[!, :distance_to_centroid] .= distances
    df[!, :distance_to_centroid_normalized] .= nd
    df[!, :centroid_channel_x] .= ctr[1]
    df[!, :centroid_channel_y] .= ctr[2]
    df[!, :centroid_channel_z] .= ctr[3]
    df[!, :centroid_object_x] .= ctrs[:, 1]
    df[!, :centroid_object_y] .= ctrs[:, 2]
    df[!, :centroid_object_z] .= ctrs[:, 3]
	# TODO compute EMST of centroids 
    return df
end

"""
    get_defaults()
    Default arguments for contact detection.
"""
function get_defaults()
    default_args = Dict()
    default_args["inpath"] = ""
    default_args["normalize"]= false
    default_args["cube-vesicle-size-ln"]=9
    default_args["cube-vesicle-sample-size"]=5
    default_args["cube-vesicle-intensity-mean"]=.2
    default_args["filtermode"]="arithmetic"
    default_args["prc"]=1.0
    default_args["denominator"]=1.0
    default_args["weighted"]=false
    default_args["sphericity"]=1.0
    default_args["nooutput"]=false
    default_args["dry-run"]=false
    default_args["noutput"]=false
    default_args["save-numerical-data"]=false
    default_args["inregex"]="*[1,2].tif"
    default_args["outpath"] = ""
    default_args["mode"] = "non-decon"
    default_args["deconvolved"] = true
    default_args["sigmas"] = "1-1-1"
    default_args["lpsigmas"] = "1-1-1"
    default_args["windowsize"] = 1
    default_args["radius"] = false
    default_args["zscore"] = 3
    default_args["volumethreshold"] = 0
    default_args["volumethresholdchannel"] = 1
    default_args["dimension"]=3
    default_args["minzslice"]=1
    default_args["alpha"]=0.05
    default_args["beta"]=0.05
    default_args["skipstats"]=false
    return default_args
end

"""
    recursive_glob(pattern, directory)
    Returns all matches files defined by `pattern` in directory, recursively
    Aimed at finding files, not directories. If your pattern matches directories, unexpected results can follow.
"""
function recursive_glob(pattern, dir)
	# Does not recurse
    files = Glob.glob(pattern, dir)
    content = readdir(dir, join=true)
    for con in content
        if isdir(con)
            files = vcat(files, recursive_glob(pattern, con))
        end
    end
    return files
end

function dtocent(ccs)
	@debug "Update to masked version"
    # cts = Images.component_centroids(ccs)
    # ac = toarray(cts[2:end])
    # m = mean(ac, dims=1)
    # ds = sqrt.(sum((ac .- m).^2, dims=2))
	ds, distances_norm, centroid, centroids = _dtocent(ccs)
    return ds
end

function _dtocent(ccs)
    cts = Images.component_centroids(ccs)
    centroids = _toarray(cts[2:end]) # Array of Nx3
    centroid = mean(centroids, dims=1)
    distances = sqrt.(sum((centroids .- centroid).^2, dims=2))
    distances_norm = _normalizemaxmin(distances)
    return distances, distances_norm, centroid, centroids
end

function _normalizemaxmin(arr)
    M = maximum(arr)
    m = minimum(arr)
    @assert all(arr .>= 0)
    @assert m <= M
    if M == m
        @debug "Max == min, degenerate normalization"
        #If the min == max, then all values are either 0 or 1. 
        if iszero(m)
            return arr # all zeros anyway
        else
            # Max == min, and it's not zero, so then normalized they're all 1
            return ones(size(arr))
        end
    else
        return (arr .- m) ./ (M-m)
    end
end


function _toarray(cts)
    if length(cts) == 0
        @error "Converting empty array"
        throw(ArgumentError("Array argument is empty"))
    end
    d = length(cts[1])
    @debug "Have $d dimensions"
    if d < 2 || d > 3
        @error "Dimensions of $d not supported, expecting 2 or 3 d."
        throw(ArgumentError("D $d not supported"))
    end
    res = zeros(Float64, length(cts), length(cts[1]))
    for (i,c) in enumerate(cts)
        if d == 3
            res[i,:] .= [c[1], c[2], c[3]]
        end
        if d == 2
            res[i,:] .= [c[1], c[2]]
        end
    end
    return res
end

function _dimg(x)
    ys = Float64.(x[:])
    if iszero(ys)
        @warn "Return NaN for zeroed image. Describing zero is unlikely what you wanted."
        return [NaN for _ in 1:8]
    end
    ys = ys[ys .> 0]
    Q1, med, Q3 = quantile(ys, [0.25, 0.5, .75])
    mx = mean(ys)
    N = length(ys)
    m2 = sum((ys .- mx).^2)/N
    m4 = sum((ys .- mx).^4)/N
    kurt = m2/m4
    m, M = minimum(ys), maximum(ys)
    return m, Q1, mx, med, Q3, M, std(ys), kurt
end


function toarray(cts)
       res = zeros(Float64, length(cts), length(cts[1]))
       for (i,c) in enumerate(cts)
           res[i,:] .= [c[1], c[2], c[3]]
       end
       return res
end

function computefeatures(fname)
	@warn "Deprecated"
    sp, qf, zs, ps = SubPrecisionContactDetection.readgld(fname)
    df_unfiltered = SubPrecisionContactDetection.reportvolumes(sp, ps, mito=sp)
    df_filtered = SubPrecisionContactDetection.reportvolumes(qf, ps, mito=sp)
    confmap = Colocalization.tomask(qf) .* (1 .- ps)
    return df_unfiltered, df_filtered, confmap
end

function readgld(fname)
    # file = jldopen(fname, "r")
    # spearman = read(file, "rawcontacts")
    # filtered = read(file, "filteredcontacts")
    # gradientcontacts = read(file, "gradientcontacts")
    # sigmap = read(file, "sigmap")
    # close(file)
	spearman = load(fname, "rawcontacts") # -> "world"
	filtered = load(fname, "filteredcontacts") # -> "world"
	gradientcontacts = load(fname, "gradientcontacts") # -> "world"
	sigmap = load(fname, "sigmap") # -> "world"
    return spearman, filtered, gradientcontacts, sigmap
end

"""
	filter_k
	Returns img where < μ + k σ is set to 0
"""
function filter_k(img, k, dropzeros=false, autotune=false, PRC=1.0, curvature=false, denominator=1.0)
	if curvature
		_img = copy(img)
		@info "Using curvature estimation"
		th = estimatecurvature(Float32.(_img[:]))
		@info "Elbow point unscaled threshold $th"
		@assert denominator > 0
		th /= denominator
		@info "Scaled threshold $th"
		_img[_img .<= th].= zero(eltype(img))
		return _img, th
	end
	return SPECHT.filter_k(img, k, dropzeros, autotune, PRC)
end

function thresholdgteq(xs, tvalue)
	ys = copy(xs)
	z = zero(eltype(xs))
	ys[xs .< tvalue].=z
	return ys
end

function getskeleton(img)
	if ! (1 < length(size(img)) < 4)
		@error "Invalid image to skeletonize"
		error("Invalid image to skeletonize.")
	end
	sk = PyCall.pyimport("skimage.morphology")
	res = sk.skeletonize_3d(Float64.(img))
	res[res .> 0] .= 1
	return Images.N0f8.(res)
end



"""
	estimatecurvature
	Estimate curvate of ys using bootstrapping.
	set seed > 0 to seed rng
"""
function estimatecurvature(ys, iters=500, samplesize=5000, seed=0)
	if seed == 0
		@warn "Not seeding RNG, non deterministic estimation"
	else
		@assert seed > 0
		Random.seed!(seed)
	end
	@info "Loading PyCall module"
    kn = PyCall.pyimport("kneed")
	@info "Sorting"
    yss = sort(ys[:])
	samplesize = min(samplesize, length(ys[:]))
    ebs = zeros(Float32, iters)
    @debug "Estimating curvature for $(size(ys)) points with $(iters) iterations and N = $(samplesize)"
    @showprogress for i in 1:iters
        yrs = Float64.(sample(yss, samplesize, replace=false, ordered=true));
        eb = kn.KneeLocator(Int32.(1:length(yrs)), Float32.(yrs), curve="convex", online=true, direction="increasing", interp_method="polynomial")
        ebs[i] = eb.elbow_y
    end
    @debug "Mean curvature $(mean(ebs)) +- $(std(ebs))"
    return mean(ebs)
end

function compute_contact_slice(target_slice, contact_slice, slice_id)
    coms = get_components(target_slice)
    NL = maximum(coms)
    indices = Images.component_indices(coms)[2:end]
    lengths = Images.component_lengths(coms)[2:end]
    contact_surfaces = zeros(size(lengths))
    for component_i in 1:NL
        # @inbounds ind = indices[component_i]
        @inbounds contact_surfaces[component_i] = sum(contact_slice[indices[component_i]])
    end
    @assert sum(target_slice) == sum(lengths)
    df = DataFrames.DataFrame(surface=lengths , surfacecontact=contact_surfaces, ratio=contact_surfaces./(lengths), componentid=1:NL, sliceid=slice_id)
    return df
    # return a dataframe
end

function getbox(boxes, image, k)
    bk = boxes[k]
    return copy(image[bk[1][1]:bk[2][1],bk[1][2]:bk[2][2],bk[1][3]:bk[2][3]])
end

"""
	imgmoment
	Returns sum(xi, yi, zi) of a 3D image
"""
function imgmoment(img)
	X, Y, Z = size(img)
	zr = zero(eltype(img))
	N = length(img[img .> zr])
	@debug "Have $N non zero voxels"
	res = zeros(Float64, N, 3)
	c = 1
	for x in 1:X, y in 1:Y, z in 1:Z
		i = img[x,y,z]
		if i > zr
			res[c,:] = [x*i, y*i, z*i]
			c += 1
		end
	end
	return res
end


"""
	shape_component
	For labeled components of an image and the i'th component, returns the eigenvalues describing its shape
"""
function shape_component(coms, image, component_index)
    boxes = Images.component_boxes(coms)[2:end]
    from, to = boxes[component_index]
    x,y,z = from
    X, Y, Z = to
    obj = copy(image[x:X, y:Y, z:Z])
    mask = coms[x:X, y:Y, z:Z]
    obj[mask .!= component_index] .= 0
    eigs = imgpca(obj)
    return eigs
end

"""
	imgpca
	Compute the principal component eigenvalues of a 3D image (centered on its non-zero mean)
	Only non-zero values are considered.
	Returns 1x3 vector of eigvals
"""
function imgpca(img)
	@assert(! all(img .== 0))
	if length(img[img .> 0]) == 1
		return [0. , 0. , 0.] #Point
	else
		rs = imgmoment(img)
		ctr = mean(rs, dims=1)
		Rc = rs .- ctr
		U, S, V = Statistics.svd(Rc)
		if size(Rc)[1] <= 2
			return S[1], S[2], 0.0 # 2xk --> 2 eigenvalues -> Line
		else
			return S
		end
		return S
	end
end


function reportvolumes(contacts, sigmap; mito=nothing)
	binarystack = contacts
    b = copy(binarystack)
    b[b .> 0] .= 1
	## Changed 3-2 connectivity
    coms = get_components(b)
    lengths = Images.component_lengths(coms)[2:end]
    indices = Images.component_indices(coms)[2:end]
    boxes = Images.component_boxes(coms)[2:end]
    N = maximum(coms)
    w=zeros(N, 15)
	@debug "Processing $N components"
	if N == 0
		@warn "NO COMPONENTS TO PROCESS"
		return nothing, nothing
	end
    @showprogress for ic in 1:N
        vals = binarystack[indices[ic]]
        valssig = sigmap[indices[ic]]
        s = sum(vals)
        gm, gs = ERGO.gmsm(vals)
        gvm, gvs = ERGO.gmsm(abs.(log.(valssig[valssig .!= 0])))
        k = kurtosis(Float64.(vals[:]))
        m, M = minimum(vals), maximum(vals)
        w[ic, 1] = s
        w[ic, 2:4] .= gm, gs, k
        w[ic, 5:6] .= m, M
        w[ic, 7:8] .= gvm, gvs
        _xy, _z, _zp = getextent(boxes[ic])
        w[ic, 9:11] .= _xy, _z, _zp
		s = slicebox(binarystack, boxes[ic])
		l1, l2, l3 = imgpca(s)
		w[ic, 13:15] .= l1, l2, l3
        @debug "Vol $s μ $gm σ $gs K $k [ $m , $M ] at extent XY $(_xy) range z $(_z) position Z $(_zp)"
    end
	w[:,12] = dtocent(coms)
	NDIST = normalizemaxmin(w[:,12])
	NZPOS = normalizemaxmin(w[:,11])
	aniso = SubPrecisionContactDetection.anisotropy.(w[:,13], w[:,14], w[:,15])
	plan = SubPrecisionContactDetection.planar.(w[:,13], w[:,14], w[:,15])
	spher = SubPrecisionContactDetection.sphericity.(w[:,13], w[:,14], w[:,15])
	## Quantify skeleton
	resadj = zeros(N, 3)
	skel = nothing
	if isnothing(mito)
		@info "No Mito -- Skipping adj info"
	else
		resadj, skel = quantify_adj_mito(mito, contacts)
	end
    df = DataFrames.DataFrame(volume=lengths, weighted=w[:,1], geometricmean=w[:,2], geometricstd=w[:,3],
    kurtosis=w[:,4], mins=w[:,5], maxs=w[:,6], meanlogsig=w[:,7], stdlogsig=w[:,8],xyspan=w[:,9],
    height=w[:,10], zposition=w[:,11], distancetocentroid=w[:,12], eig1=w[:,13], eig2=w[:,14], eig3=w[:,15],
	normalizeddistancetocentroid=NDIST, normalizedzposition=NZPOS, anisotropy=aniso, planar=plan, sphericity=spher,
	adj_mito_vol=resadj[:, 1], adj_mito_vol_fuzzy=resadj[:, 2], skeletonsurface=resadj[:,3]
	)
    return df, skel
end

function normalizemaxmin(arr)
	return (arr .- minimum(arr)) ./ (maximum(arr)-minimum(arr))
end


function sphericity(l1, l2, l3)
	if l1 == zero(eltype(l1))
        return 0
    end
    return sqrt((l2/l1) * (l3/l1))
end

function planar(l1, l2, l3)
	if l1 == zero(eltype(l1))
        return 0
    end
    return 1 - (l3/l1)
end

function anisotropy(l1, l2, l3)
    if l1 == zero(eltype(l1))
        return 0
    end
    return 1 - l2/l1
end

function slicebox(img, box)
	x, y, z = box[1]
	X, Y, Z = box[2]
	return img[x:X, y:Y, z:Z]
end

function reportvolumes2D(binarystack)
    b = copy(binarystack)
    b[b .> 0] .= 1
    coms = get_components(b)
    lengths = Images.component_lengths(coms)[2:end]
    indices = Images.component_indices(coms)[2:end]
    N = maximum(coms)
    w=zeros(N, 6)
    for ic in 1:N
        vals = binarystack[indices[ic]]
        s = sum(vals)
        gm, gs = ERGO.gmsm(vals)
        k = kurtosis(Float64.(vals[:]))
        m, M = minimum(vals), maximum(vals)
        w[ic, 1] = s
        w[ic, 2:4] .= gm, gs, k
        w[ic, 5:6] .= m, M
    end
    df = DataFrames.DataFrame(surface=lengths, weighted=w[:,1], geometricmean=w[:,2], geometricstd=w[:,3],
    kurtosis=w[:,4], mins=w[:,5], maxs=w[:,6])
    return df
end

"""
	For a bounding box, get the XY span (diagonal), Z range, and z center
"""
function getextent(box)
    xr, yr, zr = abs.(box[1] .- box[2]) .+ 1
    xy = sqrt(xr^2 + yr^2)
    return xy, zr, min(box[1][3], box[2][3]) +zr/2
end


function computesurfaces(target, contacts)
    @assert(length(size(target))==3)
    @assert(size(target)==size(contacts))
    _, _, Z = size(target)
    dfx = nothing
    # df = DataFrame(pixels=res[:,1],intensity=res[:,2], diff=res[:,3])
   # CSV.write(joinpath(OUTPATH, "test_spots_ch01_pc3ptrf.csv"), df)
    for s_i in 1:Z
        if sum(target[:,:,s_i]) ==0
            continue
        end
        @inbounds df = compute_contact_slice(target[:,:,s_i], contacts[:,:,s_i], s_i)
        if isnothing(dfx)
            dfx = df
        else
            dfx = vcat(dfx, df)
        end
    end
    return dfx
end



function process_contact_stack(tiffiles, k, w, minz=nothing, maxz=nothing; geometric=false, autotune=false, prc=1.0, curvature=false, denominator=1.0)
	## TODO 2D case, rewrite basedon 3D
    img_1 = Images.load(tiffiles[1])
    img_2 = Images.load(tiffiles[2])
    if !isnothing(minz) && !isnothing(maxz)
        @assert(1 <= minz < maxz <= size(img_1, 3))
        img_1 = img_1[:,:,minz:maxz]
        img_2 = img_2[:,:,minz:maxz]
    end
	@debug "Reuse filter_k"
    m1, s1 = mean(img_1), std(img_1)
    m2, s2 = mean(img_2), std(img_2)
    img_1f = copy(img_1)
    img_2f = copy(img_2)
    th1 = m1 + k*s1
    th2 = m2 + k*s2
	i1f, i1th = filter_k(img_1, k, geometric, autotune, prc, curvature, denominator)
	i2f, i2th = filter_k(img_1, k, geometric, autotune, prc, curvature, denominator)
	println(th1)
	println(i1th)
    img_1f[img_1 .< th1].= zero(eltype(img_1))
    img_2f[img_2 .< th2].= zero(eltype(img_2))
    sf = computecontacts(img_1, img_2, w)
    sf[(img_1f .* img_2f) .== 0].= zero(eltype(sf))
    sf[sf .> zero(eltype(sf))].= zero(eltype(sf))
    sf[sf .< zero(eltype(sf))] .= zero(eltype(sf)) .- sf[sf .< zero(eltype(sf))]
    return sf, img_1f, img_2f
end

"""
    edge_stack(image, dimensions)

    Returns the edges of image (after binarizing it), where dimensions (2, 3) indicates if the 3D dimension is to be taken into account

"""
function edge_stack(imgf, dim)
    es = copy(imgf)
    es[es .> 0].= 1
    args = (true, true, false)
    if dim == 3
        args = (true, true, true)
    end
    es = imfilter(es, Kernel.Laplacian(args))
    es[es .> 0].= 0
    es[es .< 0].= 1
    return es
end

function computecontacts(st1, st2, w)
    @assert(w > 0)
    @assert(size(st1) == size(st2))
    @assert(length(size(st1))==3)
    spear = zeros(Images.Gray{Float32}, size(st1))
    _, _, stack = size(st1)
    @threads for s in 1:stack
        @debug "Processing $s / $stack"
        _e1, _e2, _, _ = SPECHT.compute_nl(st1[:,:,s], st2[:,:,s])
        @inbounds spear[:,:,s] = sp(_e1, _e2, w)
    end
    return spear
end


function computeintensitycorrelation(st1, st2, w)
    @assert(w > 0)
    @assert(size(st1) == size(st2))
    @assert(length(size(st1))==3)
    spear = zeros(Images.Gray{Float32}, size(st1))
    _, _, stack = size(st1)
    @threads for s in 1:stack
        @debug "Processing $s / $stack"
        @inbounds spear[:,:,s] = sp(st1[:,:,s], st2[:,:,s], w)
    end
    return spear
end

function toct(img, imgl, ccs)
    _ccs = copy(img)
    _ccs[ccs .!= 0] .= 1
    _ccs[ccs .== 0] .= 0
    igl = rlap(imgl)
    return Images.colorview(Images.RGB, img, igl, _ccs)
end

function c3(vector)
    r, c = size(vector[1])
    l = length(vector)
    @assert l > 0
    result = zeros(l, r, c)
    for i in 1:l
        result[i, :, :] .= vector[i]
    end
    return result
end
"""
	expanddim(array)

	Return array of shape (size(array)..., 1)
	# Source https://stackoverflow.com/questions/42312319/how-do-i-add-a-dimension-to-an-array-opposite-of-squeeze#42312320
"""
function expanddim(array)
	S = size(array)
	return reshape(array, (S..., 1))
end

"""
	reducedim(array)

	Remove the last dimension (1) from an array
"""
function reducedim(array)
	S = size(array)
	@assert(S[end] == 1)
	reshape(array, S[1:end-1])
end

function findchannel(cs)
    channels = size(cs)[1]
    nzs = []
    @inbounds for c in 1:channels
        if any(cs[c,:,:,:] .!= 0)
            println("Found channel $c not zero")
            push!(nzs, c)
            # return c
        end
    end
    if length(nzs) == 1
        return nzs[1]
    else
        return -1
    end
end
#
# function process_stack(st1, st2, w)
#     e1 = Images.Gray{Float64}.(st1)
#     e2 = Images.Gray{Float64}.(st1)
#     e1 .= 0
#     e2 .= 0
#     spear = Images.Gray{Float64}.(st1)
#     spear .= 0.0
#     _, _, stack = size(st1)
#     _, _, _stack = size(st2)
#     @assert(stack == _stack)
#     @assert(stack >= 1)
#     # println(size(st1))
#     @threads for s in 1:stack
#         println("Processing $s / $stack")
#         _i1, _i2 = st1[:,:,s], st2[:,:,s]
#         _e1, _e2, _ep1, _ep2 = SPECHT.compute_nl(_i1, _i2) # Neg Lapl.
#         e1[:,:,s] = _ep1
#         e2[:,:,s] = _ep1
#         spear[:,:,s] = sp(_e1, _e2, w)
#     end
#     return e1, e2, spear
# end


function sp3d(img1, img2, stride; subsampleincrement=0)
    X, Y, Z = size(img1)
    spx = zeros(Float32,X,Y,Z)
	zs = zeros(Float32,X,Y,Z)
	# ts = zeros(Float32,X,Y,Z)
    @assert(size(img1)==size(img2))
    @assert(stride >= 1)
    @assert(X>stride)
    @assert(Y>stride)
    @assert(Z>stride)
	N = (X-stride*2)*(Y-stride*2)*(Z-stride*2)
	p = ProgressMeter.Progress(N, 0.5)
    @threads for x in 1+stride:X-stride
        for y in 1+stride:Y-stride
             for z in 1+stride:Z-stride
                @inbounds view1 = @view img1[x-stride:x+stride, y-stride:y+stride, z-stride:z+stride] # view saves 25% memory overhead
                @inbounds view2 = @view img2[x-stride:x+stride, y-stride:y+stride, z-stride:z+stride]
                _s, _z, _t = spearsubsample(view1, view2; increment=subsampleincrement)
                @inbounds spx[x,y,z] = _s
				@inbounds zs[x,y,z] = _z
				# @inbounds ts[x,y,z] = _t
            end
			ProgressMeter.next!(p, step=(Z-stride*2))
        end
    end
    return spx, zs, nothing
end

function compute_edges(s1, s2)
    ll1 = imfilter(s1, Kernel.Laplacian());
    ll2 = imfilter(s2, Kernel.Laplacian());
    ll1[ll1 .> 0] .= 0
    ll2[ll2 .> 0] .= 0
    return ll1, ll2
end

function sp(img1, img2, stride)
    X, Y = size(img1)
    spx = zeros(Float32,X,Y)
    @assert(size(img1)==size(img2))
    @assert(stride >= 1)
    @assert(X>stride)
    @assert(Y>stride)
    @inbounds for x in 1+stride:X-stride
        @inbounds for y in 1+stride:Y-stride
            view1 = img1[x-stride:x+stride, y-stride:y+stride]
            view2 = img2[x-stride:x+stride, y-stride:y+stride]
            srt1 = sortperm((view1[:]))
            srt2 = sortperm((view2[:]))
            # p = cor(srt1, srt2)
            spx[x,y] = cor(srt1, srt2)
        end
    end
    return spx
end

function spcor(left, right)
	@assert(size(left) == size(right))
	srt1 = sortperm(left[:])
    srt2 = sortperm(right[:])
    return cor(srt1, srt2)
end

function summarize_spots(img, imgl, concoms)
    NL = maximum(concoms)
    res = zeros(Float32, NL, 3)
    indices = Images.component_indices(concoms)[2:end]
    res[:,1]= Images.component_lengths(concoms)[2:end]
    # counts = countmap(concoms[concoms .> 0])
    @threads for component in 1:NL
        ind = indices[component]
        mskbf = sum(img[ind])
        mskbl = sum(imgl[ind])
        res[component,2] = Float32(mskbf)
        res[component,3] = Float32(mskbl)
    end
    return res
end


function binarize(img, th)
    """
    Split an array into values below, above and scaled to th.
    Returns img./th, binary mask img > th, binary mask img < 1

    """
    i = copy(img)
    @assert(! iszero(th))
    i ./= th
    et = eltype(i)
    z = zero(et)
    on = one(et)
    ones = zeros(et, size(i))
    ones[i .>= on] .= on
    ones[i .< on] .= z
    nones = zeros(et, size(i))
    nones = on .- ones
    nones[i .== z] .= z
    return i, ones, nones
end

function offset(v, low, high, tolerance)
    if v-tolerance < low
        return v + (v-tolerance) + 1
    end
    if v+tolerance > high
        return v -  ((v+tolerance)-high)
    end
    return v
end

function mcc(tp, tn, fp, fn)
    nom = (tp*tn - fp * fn)
    denom = √((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return nom/denom
end



function makespheres(radii, centers, X, Y, Z, sigma, cylaxis=nothing)
    _img1 = Images.Gray{Images.N0f16}.(zeros(X,Y,Z))
    _img2 = Images.Gray{Images.N0f16}.(zeros(X,Y,Z))
    ctr1 = omitdim(centers[1], cylaxis)
    ctr2 = omitdim(centers[2], cylaxis)
    for x in 1:X
        for y in 1:Y
            for z in 1:Z
                xs = omitdim([x y z], cylaxis)
                ed = sqrt(sum((xs - ctr1).^2))
                if ed < radii[1]
                    _img1[x,y,z] = 1.0
                else
                    _img1[x,y,z] = 0.
                end
                ed2 = sqrt(sum((xs - ctr2).^2))
                if ed2 < radii[2]
                    _img2[x,y,z] = 1.0
                else
                    _img2[x,y,z] = 0.
                end
            end
        end
    end
    _img1f = Images.imfilter(_img1, Kernel.gaussian((sigma,sigma,sigma)))
    _img2f = Images.imfilter(_img2, Kernel.gaussian((sigma,sigma,sigma)))
    return _img1, _img2, _img1f, _img2f
end

function omitdim(xs, dim=nothing)
    if isnothing(dim)
        return xs
    end
    @assert(1 <= dim <= length(xs))
    return vcat(xs[1:dim-1], xs[dim+1:end])
end


function volume_to_radius(vol)
	@assert vol > 0
    return (3/(4*π) * vol)^(1/3)
end


function radius_to_volume(radius)
	@assert radius > 0
    return (radius^3)*((4*π)/3)
end


"""
	Split a channel and its overlapping contacts based on 2 conditions.
	Used to detect faint vesicles in mito channel
	Return lower, higher, lowercontact, highercontact (masks)
	treshold : (weighted) volume
	intensityfilter : apply a weak filter to split faintly connected parts
	sphericity_threshold : 1 = spherical
"""
function filter_channels(channel, interact, threshold; weighted=false, intensityfilter=false, sphericity_threshold=-1)
	img = copy(channel)
	contacts = copy(interact)
	sphericity_threshold != -1 ? @assert(false) : @debug "Skipping sphr test"
	if intensityfilter
		@error "Using intensity filter --> UNSTABLE"
		img = filterintensity(img)
	end
	lower = Colocalization.aszero(img)
	higher = copy(lower)
	IM = Colocalization.tomask(img)
	channel_components = get_components(IM)
	lengths = Images.component_lengths(channel_components)[2:end]
	@showprogress for (ith, component_indices) ∈ enumerate(Images.component_indices(channel_components)[2:end])
		component = img[component_indices]
		V = weighted ? sum(component) : lengths[ith]
		if (V <= threshold)
			lower[component_indices] .= 1
		else
			higher[component_indices] .= 1
		end
	end
	@assert sum(lower .* higher) == 0
	lower_contacts = Colocalization.aszero(contacts)
	higher_contacts = Colocalization.aszero(contacts)
	CM = Colocalization.tomask(contacts)
	contacts_components = get_components(CM)
	@showprogress for (jth, contact_indices) ∈ enumerate(Images.component_indices(contacts_components)[2:end])
		if any(lower[contact_indices] .> 0)  # Touching vesicle --> v contact
			lower_contacts[contact_indices] .= 1
		else ## otherwise non v contact
			higher_contacts[contact_indices] .= 1
		end
	end
	@assert sum(lower_contacts) + sum(higher_contacts) == sum(CM)
	return lower, higher, lower_contacts, higher_contacts
end



function computesphericity(mask)
	coms = get_components(mask)
	S = zeros(maximum(coms))
	for (i,(c,b)) in enumerate(zip(Images.component_indices(coms)[2:end], Images.component_boxes(coms)[2:end]))
		s = slicebox(mask, b)
		S[i] = sphericity(imgpca(s)...)
	end
	return S
end


function filterintensity(intensity, d=4)
	th = estimatecurvature(SPECHT.nz(Float64.(intensity[:])), 500, 5000, 42)
	filtered = copy(intensity)
	# @error "Threshold $(th/d)"
	filtered[filtered .< th/d].=0
	return filtered
end

function normnz(img)
	cimg = copy((img))
	min = minimum(SPECHT.nz(cimg))
	max = maximum(SPECHT.nz(cimg))
	cimg[cimg .> 0] .= (cimg[cimg .> 0] .- min)./(max-min)
	return cimg
end


function masktoindices3(A)
	X, Y, Z = size(A)
	return [[x, y, z] for z in 1:Z, y in 1:Y, x in 1:X if A[x,y,z]>zero(eltype(A))]
end
function masktoindices2(A)
	X, Y = size(A)
	return [[x, y] for y in 1:Y, x in 1:X if A[x,y,]>zero(eltype(A))]
end

function masktoindices(A)
	nd = length(size(A))
	if nd == 2
		return masktoindices2(A)
	end
	if nd == 3
		return masktoindices3(A)
	end
	@assert(false)
end

end # module
