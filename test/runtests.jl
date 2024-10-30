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
import PyCall
using SubPrecisionContactDetection
import ERGO
using Test
import Colocalization
using Images, Colors, DataFrames, CSV, Statistics, LinearAlgebra #Imageview crashes in headless
import Glob
import ImageMagick
import ImageFiltering
import Random
import Match
using Distributions

@testset "SubPrecisionContactDetection.jl" begin

    # HTests updated their API, which now breaks. I'm not using this function, so disable until I have time to fix it.
	@testset "gdf" begin
        @warn "Fix HTests"
	# 	m, M = getci([100,101,100], 0.05)
	# 	@test m < M
	end

	@testset "fc" begin
		A = zeros(3, 100, 100, 100)
		@test findchannel(A) == -1
		A[2,:,:,:] .= 1
		@test findchannel(A) == 2
	end

    @testset "3dshape" begin
        Random.seed!(42)
        i = zeros(10, 10, 10)
        i[5:6, 5:6, 5:6] .= rand()
        i[1:1, 1:1, 1:1] .= rand()
        coms = Images.label_components(i)
        eigs = shape_component(coms, i, 1)
        @test length(eigs) == 3
        @test eigs[1] >= eigs[2] >= eigs[3]
        eigs = shape_component(coms, i, 2)
        @test length(eigs) == 3
        @test eigs[1] >= eigs[2] >= eigs[3]
        @test eigs[1] > 0
        @test isapprox(eigs[1], eigs[2], atol=1e-3)
        @test isapprox(eigs[2], eigs[3], atol=1e-3)
        i = zeros(10, 10, 10)
        i[5:5, 5:8, 5:5] .= rand()
        coms = Images.label_components(i)
        eigs = shape_component(coms, i, 1)
        @test isapprox(eigs[2], 0, atol=1e-1)
        @test isapprox(eigs[1], 1, atol=0.2)
    end

	@testset "s2" begin
		A = zeros(100, 100)
		r=sp2d(A, A, 1)[1]
		unique(r) == [0, 1]
	end

    @testset "recursiveglob" begin
        testdir = "testglob"
        testfile = "1.tif"
        testpattern = "*1.tif"
        if isdir(testdir)
            rm(testdir, recursive=true)
        end
		mkdir(testdir)
        mkdir(joinpath(testdir, "1"))
        mkdir(joinpath(testdir, "1", "2"))
        touch(joinpath(testdir, testfile))
        touch(joinpath(testdir, "1", testfile))
        touch(joinpath(testdir, "1", "2", testfile))
        fs = Glob.glob(testpattern, testdir)
        fx = recursive_glob(testpattern, testdir)
        @test length(fs) == 1
        @test length(fx) == 3
        for f in fs
            @test isfile(f)
        end
        if isdir(testdir)
            rm(testdir, recursive=true)
        end
	end

    @testset "comb" begin
        t = mktempdir()
		Images.save(joinpath(t, "a01.tif"), rand(200, 200))
        Images.save(joinpath(t, "b02.tif"), rand(200, 200))
        Images.save(joinpath(t, "c03.tif"), rand(200, 200))
        fs = Glob.glob("*.tif", t)
        @test length(fs) == 3
        es = endings(fs)
        @test es[1] == 1
        @test es[2] == 2
        @test es[3] == 3
        cs, cis = combines(es)
        @test cs[1] == [1,2]
        @test cs[2] == [1,3]
        @test cs[3] == [2,3]
        x, xis = combines([1,2,3])
        @test x == xis
		rm(t, recursive=true) 
    end

    @testset "mulitchannel" begin
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
        rx = "*[1,2,3].tif"
        paths, regexes = buildregex(subd, rx)
        @test length(paths) == 3
        @test length(regexes) == 3
        @test paths == ["1--2", "1--3", "2--3"]
        @test regexes == ["*[1,2].tif", "*[1,3].tif", "*[2,3].tif"]
		rm(t, recursive=true) 
    end

    @testset "meta" begin    
        args = get_defaults()
        save_meta("test.json", args)
        fargs = read_meta("test.json")
        @test args == fargs
    end

    @testset "def" begin
        defs = get_defaults()
        @test length(defs) == 30
    end

    @testset "2channelrefactor" begin
        # defs = get_defaults()
        t = mktempdir()
        c1 = Images.N0f16.(zeros(100, 100, 100))
        c1[50:60, 50:60, 50:60] .= 1
        c2 = Images.N0f16.(zeros(100, 100, 100))
        c2[58:70, 58:70, 58:70] .= 1
        sigma = 1.0
        c1f = Images.imfilter(c1, Kernel.gaussian((sigma,sigma,sigma)))
        c2f = Images.imfilter(c2, Kernel.gaussian((sigma,sigma,sigma)))
        clamp01!(c1f)
        clamp01!(c2f)
        Images.save(joinpath(t, "1.tif"), c1f)
        Images.save(joinpath(t, "2.tif"), c2f)
        defs = get_defaults()
        defs["zscore"] = 0.5
        defs["inpath"] = t
        q = mktempdir()
        defs["outpath"] = q
        two_channel_contacts(defs)
        @test length(Glob.glob("*.tif", q)) == 6
        rm(t, recursive=true)
        rm(q, recursive=true)
    end

	@testset "iq" begin

		A = zeros(100, 100)
		r=sp2d(A, A, 1)[1]
		unique(r) == [0, 1]
		using Images
		t = mktempdir()
		Images.save(joinpath(t, "x.tif"), rand(200, 200))
		r=reportimagequality(joinpath(t, "x.tif"))
		all(r .> 0)
		rm(t, recursive=true)

	end

    @testset "ce" begin
		A = zeros(100, 100)
		l1, l2 = compute_edges(A, A)
		@test l1 == l2
	end

	@testset "api" begin
		Random.seed!(42)
		A = rand(100, 100, 20)
		B = rand(100, 100, 20)
		t=mktempdir()
		Images.save(joinpath(t, "1.tif"), A)
		Images.save(joinpath(t, "2.tif"), B)
		tiffiles = [joinpath(t, "$x.tif") for x in 1:2]
		sf, img_1f, img_2f, zs, _, ll1, ll2, qf, rawsf=process_contact_stack3d(tiffiles, 3, 2)
		@test sum(sf) < 1
	end

	@testset "normlin" begin
		a = [0 0; 0.5 1]
		b = normalize_linear(a)
		@test maximum(b) == 1
		@test minimum(b) == 0
		@test sum(b) == 1
		Random.seed!(42)
		for _ in 1:100
			a = rand(2,2)
			b = normalize_linear(a)
			@test minimum(b) == 0
			@test maximum(b) == 1
		end
	end

	@testset "RDSTACK" begin
		using Random
		Random.seed!(42)
		xs = rand(10, 10)
		ys = expandstack(xs, 2)
		zs = reducestack(ys, 2)
		@test all(zs .== xs)
		@test all(ys[:,:,1] .== xs)
	end

	@testset "edgefromborder" begin
		using Random
		Random.seed!(42)
		for i in 1:20
			xs = zeros(30, 30, 30)
			xs[11:20, 11:20, 11:20] .= rand(10, 10, 10)
			sx = copy(xs)
			r = edge_from_border(xs, 1)
			q = edge_from_border(xs, 2)
			@test sum(q) > sum(r)
			@test all(xs .== sx)
		end
	end


	@testset "normimg" begin
		using Random
		Random.seed!(42)
		ref = rand(10, 10, 10)
		for _ in 1:10
			c2 = rand(10, 10, 10) ./ 3
			_, i2 = normalize_channel(ref, c2)
			@test sum(c2) < sum(i2)
		end
	end

	@testset "filtermito" begin
		using Random
		Random.seed!(42)
		a = zeros(10, 10, 10)
		a[3:6, 3:6, 3:6] .= 0.5
		filtered, mint, logsize = filtermito(a, a, 100, 100)
		@test size(logsize)[1] == 1
		@test sum(filtered) == 0
		@test mint[1] == 0.5
		@test 4 < logsize[1] < 5
		filtered, mint, logsize = filtermito(a, a, 0, 0)
		@test size(logsize)[1] == 1
		@test sum(filtered) == 32
		@test mint[1] == 0.5
		@test 4 < logsize[1] < 5
	end

	@testset "walk_cube" begin
		a = zeros(10, 10, 10)
		a[3:6, 3:6, 3:6] .= 0.5
		df, b1, b2 = walk_cube(a, a, a, 5, 5, 5)
		@test sum(df[!, :mitvol]) == sum(a)*2
		@test sum(df[!, :mitosum]) == sum(a)
		@test size(df)[1] == 8
		 # ["mitvol", "mitosum", "contactvol", "contactsum", "ncontacts", "mtsurface", "ctsurface"]
		@test all(df[!, :mtsurface] .== df[!, :ctsurface])
		@test all(sum(df[!, :mtsurface]) == sum(Colocalization.tomask(b1)))
		@test all(sum(df[!, :ctsurface]) == sum(Colocalization.tomask(b2)))
		df, b1, b2 = walk_cube(a, Colocalization.aszero(a), a, 5, 5, 5)
		@test sum(b2) == 0
		@test sum(b1) > sum(b2)
		@test sum(df[!, :ctsurface]) == 0
		@test !all(df[!, :mtsurface] .== df[!, :ctsurface])
		Random.seed!(42)
		for _ in 1:100
			a = abs.(rand(10, 10, 10) .- 0.5)
			df, b1, b2 = walk_cube(a, a, a, 5, 5, 5)
			@test sum(df[!, :mitvol]) == sum(Colocalization.tomask(a))
			@test isapprox(sum(df[!, :mitosum]), sum(a), atol=0.001)
			@test size(df)[1] == 8
			@test all(sum(df[!, :mtsurface]) == sum(Colocalization.tomask(b1)))
			@test all(sum(df[!, :ctsurface]) == sum(Colocalization.tomask(b2)))
		end
		for _ in 1:100
			a = abs.(rand(10, 10, 10) .- 0.5)
			df, b1, b2 = walk_cube(a, ERGO.aszero(a), a, 5, 5, 5)
			@test iszero(b1)
			@test sum(df[!, :mitvol]) == sum(Colocalization.tomask(a))
			@test isapprox(sum(df[!, :mitosum]), sum(a), atol=0.001)
			@test size(df)[1] == 8
		end
		a = zeros(10, 10, 10)
		a[3:6, 3:6, 3:6] .= 0.5
		bd = edge_from_border(a, 1)
		c = copy(bd)
		ct = clipr(c, 0, 0.2)
		@test sum(ct) < sum(bd)
		df, b1, b2 = walk_cube(a, ct, ct, 5, 5, 5)
		@test sum(df[!,:mtsurface]) > sum(df[!,:ctsurface])
	end

	@testset "PS" begin
		st = "1-2-3"
		sig = parsesigmas(st)
		@test sig == [1, 2, 3]
		st = "1.0-2.0-3.12"
		sig = parsesigmas(st)
		@test sig == [1.0, 2.0, 3.12]
	end

	@testset "dims" begin
		using Random
		Random.seed!(42)
		xs = rand(10, 10, 10)
		ys = xs[:,:,5:5]
		yss = reducedim(ys)
		# 10, 10
		@test length(size(yss)) == 2
		yps = expanddim(yss)
		# 10, 10, 1
		@test length(size(yps)) == 3
		@test all(yps[:] .== ys[:])
	end

	@testset "glew" begin
		a = zeros(10, 10, 10)
		g = gerode(a)
		@test iszero(g)
		a[5,5,5] = 1
		g = gerode(a)
		@test iszero(g)
		g = dropleq(a, 2)
		@test iszero(g)
		g = dropleq(a, 0)
		@test sum(a) == 1
	end

	@testset "minr" begin
		ζ = 0.01
		N = 125
		minr = compute_min_r_for_sample_corr(N, ζ, ζ)
		@test isapprox(minr, 0.435, atol=0.01)
		ζ = 0.05
		minr2 = compute_min_r_for_sample_corr(N, ζ, ζ)
		@test isapprox(minr2, 0.340, atol=0.01)
		@test minr2 < minr
		@test minr > 0
		minr3 = compute_min_r_for_sample_corr(60, ζ, ζ)
		@test minr2 < minr3
	end

	@testset "sphr" begin
		mk = zeros(10, 10, 10)
		mk[4:5, 4:5, 4:5].=1
		S = computesphericity(mk)
		@test length(S) == 1
		@test all(S .> 0.5)
	end

    @testset "spc" begin
        Random.seed!(42)
        RG = colorview(RGB, rand(10,10,10), zeros(10,10,10), zeros(10,10,10))
        R = splitchannel(RG)
        @test ! iszero(R)
        RG = colorview(RGB, zeros(10,10,10), rand(10,10,10), zeros(10,10,10))
        R = splitchannel(RG)
        @test ! iszero(R)
        RG = colorview(RGB, zeros(10,10,10), zeros(10,10,10), rand(10,10,10))
        R = splitchannel(RG)
        @test ! iszero(R)
        RG = colorview(RGB, zeros(10,10,10), zeros(10,10,10), rand(10,10,10))
        QR = copy(RG)
        R = splitchannel(RG)
        @test QR == RG
        RZ = rand(10, 10, 10)
        Z = splitchannel(RZ)
        @test Z == RZ
    end

    @testset "filterchannels" begin
		a = zeros(10,10,10)
		a[5:7,5:7,5:7] .= 1
		c = zeros(10,10,10)
		c[5:7,5:7,5:7] .= 1
		L, H, CL, CH = filter_channels(a, c, 26, weighted=false)
		@test sum(L) + sum(H) == sum(a)
		@test sum(CL) + sum(CH) == sum(c)
		@test iszero(sum(L))
		@test iszero(sum(CL))
		@test !iszero(sum(H))
		@test !iszero(sum(CH))
		L, H, CL, CH = filter_channels(a, c, 28, weighted=false)
		@test sum(L) + sum(H) == sum(a)
		@test sum(CL) + sum(CH) == sum(c)
		@test !iszero(sum(L))
		@test iszero(sum(H))
		Random.seed!(42)
		for _ in 1:100
			channel = rand(100, 100, 100)
			channel[channel .> 0.6] .= 0
			c = rand(100, 100, 100)
			channel[channel .> 0.7] .= 0
			L, H, CL, CH = filter_channels(channel, c, Int(round(rand()*10)), weighted=false)
			@test sum(L) + sum(H) == sum(Colocalization.tomask(channel))
		end
    end

    @testset "fmi" begin
        Random.seed!(42)
        for _ in 1:10
            a = rand(100,100)
            ia = iterativemedian(a, 1, 1)
            @test snr(a) <= snr(ia)
        end
    end

	@testset "skel" begin
		A = zeros(5, 5, 5)
		A[2:4, 2:4, 2:4] .= 1
		SA = getskeleton(A)
		@test sum(A) > 0
		A = zeros(5, 5, 5)
		A[3:5, 3, 3] .= 1
		SA = getskeleton(A)
		@test sum(SA)==3
	end

	@testset "qam" begin
		A = zeros(5, 5, 5)
		A[2:4, 2:4, 2:4] .= 1
		SK = getskeleton(A)
		res, skel =  quantify_adj_mito(A, A)
		@test !isnothing(res)
		@test all(skel .== SK)
		@test size(res)[1] == 1
		@test res[1,1] == sum(A)
		B = zeros(10, 10)
		res, skel =  quantify_adj_mito(B, B)
		@test isnothing(res)
		A = zeros(5, 5, 5)
		B = zeros(5, 5, 5)
		A[2:4, 2:4, 2:4] .= 1
		B[3:5, 3:5, 3:5] .= 0.5
		SK = getskeleton(A)
		res, skel =  quantify_adj_mito(A, B)
		@test all(skel .== skel)
		A = zeros(5, 5, 5)
		B = zeros(5, 5, 5)
		A[2:4, 2:4, 2:4] .= 1
		B[3:5, 3:5, 3:5] .= 0.5
		SK = getskeleton(A)

		res, skel =  quantify_adj_mito(zeros(5,5,5), A)
		@test res == [0.0 0.0 3.0]
		A = zeros(5, 5, 5)
		B = zeros(5, 5, 5)
		A[2:4, 2:4, 2:4] .= 1
		B[3:5, 3:5, 3:5] .= 0.5
		SK = getskeleton(A)
		res, skel =  quantify_adj_mito(B, A)
		@test res == [27.0 13.5 3.0]
	end

    @testset "radiusvolume" begin
        Random.seed!(42)
		for _ in 1:1000
			r = rand()*10
			@test isapprox(volume_to_radius(radius_to_volume(r)), r)
		end
		vp = 129.87878804533656
		@test isapprox(vp, radius_to_volume(π))
		@test isapprox(π, volume_to_radius(radius_to_volume(π)))
    end

    @testset "ab" begin
        alpha = 0.05
        beta = 0.25
        Random.seed!(42)
        for _ in 1:10
            a = rand(100,100)
            ps = rand(100, 100)
            b = alphabeta!(copy(a), ps, alpha, beta)
            @test sum(b) < sum(a)
            @test all(b[b .> 0] .>= beta)
        end
    end

    @testset "spcor" begin
        @test spcor([1 1], [2 2])==1
        @test spcor([1 1], [1 1])==1
        @test spcor([1 2], [2 1])==-1
    end

    @testset "denoise" begin
        img = zeros(3,3,3)
        denoise(img, 2)
        img = zeros(3,3)
        denoise(img, 2)
        sigmas = [1, 1, 1]
        for option in ["harmonic", "geometric", "median", "gaussian"]
            img = zeros(3,3,3)
            nimg=denoise(img, 2, option, sigmas)
            @test ! any(isnan.(nimg))
            img = rand(10,10,10)
            cimg = copy(img)
            nimg = denoise(img, 2, option, sigmas)
            @test sum(nimg) != sum(cimg)
            @test all(cimg .== img)
            @test ! any(isnan.(nimg))
        end
    end

    @testset "festimparam" begin

        xs = [1 2 3 4]
        ys = [2 3 4 5.1]
        testx = [0.9 1 1.5 2 2.5 3 3.5 4 4.1]
        expx =  [2 2 2.5 3 3.5 4 4.55 5.1 5.1]

        for (tx, ex) in zip(testx, expx)
        	@debug "$(tx) --> $(ex) ? "
        	yest = festimparam(xs, ys, tx)
        	@debug " $(yest)"
        	@test yest == ex
        end
    end

    @testset "snr" begin
        A = [1 2 3]
        @test snr(A) == 2.0
        @test isnan(snr([0 0]))
        @test isinf(snr([1 1]))
    end

    @testset "filter" begin
        TG = zeros(4, 4, 4)
        FT = zeros(4, 4, 4)
        TG[2,2,2]=1
        FT[2,2,2]=1
        FT[3,2,2]=1
        TG[1,1,1]=1
        FP, FN, TP, TN, fp, fn, tp, tn = SubPrecisionContactDetection.ratefilter(TG, FT)
        @test sum(FP) == 1
        @test sum(TP) == 1
        @test sum(TN) == 61
        @test sum(FN) == 1
        @test fp == tp
        @test tn == 1-fn
    end

    @testset "curve" begin
        d = Distributions.Normal()
        Random.seed!(42)
        ys = rand(d, 1000)
        thx = SubPrecisionContactDetection.estimatecurvature(ys, 10, 100, 0)
        thy = SubPrecisionContactDetection.estimatecurvature(ys, 10, 100, 0)
        @test thx != thy
        thx = SubPrecisionContactDetection.estimatecurvature(ys, 10, 100, 1)
        thy = SubPrecisionContactDetection.estimatecurvature(ys, 10, 100, 1)
        @test thx == thy
    end

    # function thresholdgteq(xs, tvalue)
    @testset "threshold" begin
        d = Distributions.Normal()
        Random.seed!(42)
        xs = rand(d, 20, 20)
        mu = Statistics.mean(xs)
        ys = SubPrecisionContactDetection.thresholdgteq(xs, mu)
        @test all(ys[xs .< mu] .== zero(eltype(xs)))
    end


    @testset "clampt" begin
        xs = -.3:.01:.3 |> collect
        ys = smoothclamp(xs)
        ysk = smoothclamp(xs, 3)
        @test all(0.0 .<= ys .<= 1.0)
        @test all(0.0 .<= ysk .<= 1.0)
    end

    @testset "clampt2" begin
        a = [-1 -.5 -.4 0 .2 .5]
        th = .3
        resc = clampt(a, th)
        @test all(0.0 .<= resc .<= 1.0)
        @test resc[1] == 0
        @test resc[4] == 0.5
        @test resc[6] == 1
    end

    @testset "DTD" begin
        diag_matrix = [1 0 0; 0 4 0; 0 0 1]
        eig = eigvals(diag_matrix)
        eigv = eigvecs(diag_matrix)
        λ = zeros(1,1,1,3)
        V = zeros(1,1,1,3,3)
        λ[1,1,1,:] = eig
        V[1,1,1,:, :] = eigv
        field, mg = dtd_to_field(λ, V)
        @test field[1,1,1,:] == [0.0, 1.0, 0.0]
        @test mg[1,1,1] == 4.0
    end

    @testset "disttocentroid" begin
        A = zeros(4,4,4)
        A[1,1,1] = 1
        A[4,4,4] = 1
        cc = Images.label_components(A)
        a, b = dtocent(cc)
        @test a == b
        @test isapprox(a, 2.598076211353316)
    end

    @testset "imgpca" begin
        A = zeros(10,10,10)
        SX = sphere3d([3,3,3], 2, 10,10,10, 0)
        SubPrecisionContactDetection.indicestomask!(SX, A)
        λ1, λ2, λ3 = SubPrecisionContactDetection.imgpca(A)
        @test λ1 >= λ2 >= λ3
        @test isapprox(λ1, λ2)
        @test isapprox(λ1, λ3)
        A = ones(2, 1, 1)
        λ1, λ2, λ3 = SubPrecisionContactDetection.imgpca(A)
        @test λ1 >= λ2 >= λ3
        @test isapprox(λ1, 0.707, atol=0.1)
        @test isapprox(λ2, 0)
        @test isapprox(λ3, 0)
        A = ones(1, 1, 1)
        λ1, λ2, λ3 = SubPrecisionContactDetection.imgpca(A)
        @test λ1 >= λ2 >= λ3
        @test isapprox(λ1, 0)
        @test isapprox(λ2, 0)
        @test isapprox(λ3, 0)
        A = zeros(3,3,3)
        A[2, 2, 2] = 1
        λ1, λ2, λ3 = SubPrecisionContactDetection.imgpca(A)
        @test λ1 >= λ2 >= λ3
        @test isapprox(λ1, 0)
        @test isapprox(λ2, 0)
        @test isapprox(λ3, 0)
    end

    @testset "sigtest" begin
        xs = -3:.1:3 |> collect
        p, c, s = z_to_significance(xs)
        @test length(s[s.>0]) == 14
        @test all(c[30:end] .> p[30:end])
        @test isapprox(maximum(p), (sqrt(2*π))^-1 *  ℯ^-(0))
        @test isapprox(maximum(c), 1, atol=0.002)
        _pd, _cd, _msk = z_to_significance([1,2,3,Inf])
        @test _msk[4] == 0
        @test _pd[4] == 0
        @test _cd[4] == 1
    end

    @testset "maggrad" begin
        ll2 = Random.rand(100,100,100)
        δ2x, δ2y, δ2z = Images.imgradients(ll2, KernelFactors.ando3);
        mg2 = sqrt.(δ2x.^2 .+ δ2y.^2 .+ δ2z.^2)
        mg2p = gradientmagnitude([δ2x, δ2y, δ2z])
        mg23p = magnitudegradients(δ2x, δ2y, δ2z)
        @test isapprox(mg2,mg2p)
        @test isapprox(mg2,mg23p)
        ll2 = Random.rand(100,100)
        δ2x, δ2y = Images.imgradients(ll2, KernelFactors.ando3);
        mg2 = sqrt.(δ2x.^2 .+ δ2y.^2)
        mg2p = gradientmagnitude([δ2x, δ2y])
        @test isapprox(mg2,mg2p)
    end

    @testset "SCM" begin
        A = zeros(3, 3)
        A[1, 1] = 1
        A[3, 3] = 1
        r  = SubPrecisionContactDetection.scoremasks(A, A, A, A, A, A)
        @test r["1"] == [1, 1, 2, 2]
        @test length(keys(r)) == 7
    end


    @testset "spear" begin
        A = ones(3,3,3)
        B = copy(A)
        B[2] = 4
        @test spear(A, A)[1] == 1.0
        @test spear(A, B)[1] <= 1.0
        @test spearsubsample(A, B, increment=1)[1] <= 1.0
        @test spearsubsample(A, A, increment=1)[1] == 1
    end

    @testset "mtindex" begin
        A = zeros(5,5,5)
        B = zeros(2,2)
        A[1,1,1] = 1
        B[1,1] = 1
        @test masktoindices(A) ==[[1,1,1]]
        @test masktoindices(B) ==[[1,1]]
        inds = [[1,2], [2,2]]
        arr = zeros(Float64, 5, 5)
        @test sum(indicestomask!(inds, arr)) == 2.0
    end

    @testset "c3" begin
        a, b, c = (3, 4, 5)
        A = zeros(a, b, c)
        A[1,2,3] = 42
        V =  []
        for i in 1:a
            push!(V, A[i, :, :])
        end
        _c3 = SubPrecisionContactDetection.c3(V)
        @test isapprox(_c3[1,2,3], 42)
        @test size(_c3) == (a, b, c)
    end

    @testset "boxk" begin
        A = zeros(10, 10, 10)
        A[2:3, 4:5, 5:10] .= 1
        coms = Images.label_components(A)
        lengths = Images.component_lengths(coms)[2:end]
        indices = Images.component_indices(coms)[2:end]
        boxes = Images.component_boxes(coms)[2:end]
        k = 1
        boxk = getbox(boxes, A, k)
        @test size(boxk) == (2, 2, 6)
    end


    @testset "vol" begin
        a = Gray{N0f16}.(zeros(100, 100, 100))
        a[2:4, 2:4, 2:4] .= .5
        sm = Random.rand(size(a)...)/1000
        df, sk = reportvolumes(a, sm; mito=a)
        @test size(df) == (1,24)
        @test round(df.eig1[1]) <= df.eig2[1]
        @test round(df.eig1[1]) <= df.eig3[1]
        @test isapprox(df.eig1[1], df.eig3[1])
        @test isnan(df.normalizeddistancetocentroid[1])
        @test df.distancetocentroid[1]==0
        qdf = reportvolumes2D(a)
        @test size(qdf) == (1,7)
        @test df.volume[1] == 27
        @test isapprox(df.weighted[1], 27/2, atol=.01)
        @test isapprox(df.geometricmean[1], .5, atol=.01)
        @test isapprox(df.geometricstd[1], 1.0, atol=.01)
        @test isnan(df.kurtosis[1])
        @test isapprox(max(df.maxs[1]), .5, atol=.01)
        @test isapprox(min(df.mins[1]), .5, atol=.01)
        @test isapprox(df.distancetocentroid[1], 0)
        a = Random.zeros(10,10,10)
        b = Random.ones(10,10,10)
        res, _ = reportvolumes(a, b; mito=a)
        @test isnothing(res)
    end

    @testset "volext" begin
        a = Gray{N0f16}.(zeros(100, 100, 100))
        a[2:4, 2:4, 2:4] .= .5
        a[80, 80, 80] = .5
        a[20:22, 20:23, 20:26] .= .5
        sm = Random.rand(size(a)...)/1000
        df, _ = reportvolumes(a, sm, mito=a)
        @test isapprox(df.distancetocentroid[1], 55.3315, atol=0.01)
        @test size(df) == (3,24)
        @test round(df.eig1[1]) <= df.eig2[1]
        @test round(df.eig1[1]) <= df.eig3[1]
        @test isapprox(df.eig1[1], df.eig3[1])
        @test df.volume[1] == 27
        @test df.volume[2] == 3*4*7
        @test isapprox(df.eig1[2], 9, atol=.5)
        @test isapprox(df.eig2[2], 5, atol=.5)
        @test isapprox(df.eig3[2], 3, atol=1)
    end

    @testset "normmaxmin" begin
        A = Random.rand(10) * 200
        An = SubPrecisionContactDetection.normalizemaxmin(A)
        @test all(An .<= 1)
        @test all(An .>= 0)
    end

    @testset "featurespca" begin
        l1, l2, l3 = 10, 10, 10
        @test isapprox(sphericity(l1, l2, l3), 1)
        @test isapprox(planar(l1, l2, 0), 1)
        @test isapprox(planar(l1, l2, l3), 0)
        @test isapprox(anisotropy(l1, 0, 0), 1)
        l1, l2, l3 = 0, 0, 0
        @test isapprox(planar(l1, l2, l3), 0)
        @test isapprox(sphericity(l1, l2, l3), 0)
        @test isapprox(anisotropy(l1, l2, l2), 0)
    end

    @testset "fk" begin
        img = rand(100, 100, 10)
        u, s = Statistics.mean(img), Statistics.std(img)
        for k in 1.:4.
            ik, th = filter_k(img, k)
            @test all(ik[ik .> 0] .>= th)
            @test isapprox(u+k*s, th)
        end
        img = rand(100)
        img[25:45] .= 0.0
        ik, th = filter_k(img, 1, false)
        ikz, thz = filter_k(img, 1, true)
        @test thz > th
        @test sum(ik) > sum(ikz)
        img = 1.0 .+ rand(100)
        ik, th = filter_k(img, 1, false)
        ikz, thz = filter_k(img, 1, true)
        @test isapprox(thz, th, atol=0.02)
        sikz, sthz = filter_k(img, 1, true, true, 1.0)
        s2ikz, s2thz = filter_k(img, 1, true, true, 2.0)
        sikz, sthz = filter_k(img, 1, true, true, 1.0, true, 1.0)
        s2ikz, s2thz = filter_k(img, 1, true, true, 2.0, true, 2.0)
        @test sthz > s2thz
        @test sum(sikz) < sum(s2ikz)
    end

    @testset "edgestack" begin
        _data = zeros(10,10,10)
        _data[5:7,5:7,5:7] .= rand()
        es = edge_stack(_data, 3)
        @test all(es[5,5,5:7] .== 1)
        @test all(es[6,6,5] .== 1)
        @test all(es[6,6,7] .== 1)
        @test all(es[6,6,6] .== 0)
        es = edge_stack(_data, 2)
        @test all(es[5,5,5:7] .== 1)
        @test all(es[6,6,5] .== 0)
        @test all(es[6,6,7] .== 0)
        @test all(es[6,6,6] .== 0)
    end

    @testset "contact" begin
        _data = zeros(10,10,3)
        _data[5:7,5:7,:] .= 1
        contact = zeros(10, 10, 3)
        contact[5:7,5,2] .= 1
        df = compute_contact_slice(_data[:,:,1], contact[:,:,1], 1)
        @test all(df.ratio .== 0)
        df = compute_contact_slice(_data[:,:,2], contact[:,:,2], 2)
        @test sum(df.ratio)==1/3
        dfx = computesurfaces(_data, contact)
        @test dfx.ratio[1] == 0
        @test dfx.ratio[2] == 1/3
        @test dfx.ratio[3] == 0
    end

    @testset "correlation" begin
        X, Y, Z = 100, 100, 10
        Random.seed!(42)
        img1 = rand(X,Y,Z)
        img2 = rand(X,Y,Z)
        for w in [1,2,3]
            si = computeintensitycorrelation(img1, img2, w)
            sc = computecontacts(img1, img2, w)
            @test all(-1 .<= si .<= 1)
            @test all(-1 .<= sc .<= 1)
        end
    end

    @testset "sp" begin
        X, Y = 10, 10
        img1 = zeros(X,Y)
        stride = 2
        for x in 1:X
            for y in 1:Y
                img1[x,y] = x+y
            end
        end
        rsp = sp(img1, img1, stride)
        @test rsp[5,5] == 1
        rsp = sp(img1, 0 .- img1, stride)
        @test all(rsp .<= 0)
        @test rsp[1+stride-1,1+stride-1] == 0

    end

    # @testset "GLD" begin
    #     Random.seed!(42)
    #     for _ in 1:100
    #         a = Random.rand(20, 20, 20)
    #         a[a.<.5] .= 0
    #         b = Random.rand(20, 20, 20)
    #         b[b.<.6] .= 0
    #         c = Random.rand(20, 20, 20)
    #         d = Random.rand(20, 20, 20)
    #         savegld(joinpath(tempdir(),"test.jld2"), a, b, c, d)
    #         A, B, C, D = readgld(joinpath(tempdir(),"test.jld2"))
    #         @test eltype(A) == Float64
    #         @test A == a
    #         @test B == b
    #         @test C == c
    #         @test D == d
	# 		savegld(joinpath(tempdir(),"test.jld2"), a, b, c, d; tp=Images.N0f8)
    #         A, B, C, D = readgld(joinpath(tempdir(),"test.jld2"))
    #         @test eltype(A) == Images.N0f8
    #     end
    # end

    @testset "RN" begin
        A = [1 0 -1 0]
        am = retainnegative(A)
        @test am[3] == 1
        @test all(am[A .>= 0] .== 0)
    end


    @testset "getextent" begin
        A = zeros(5, 5, 5)
        A[2,2:4,2:3] .= Random.rand(3,2)
        binarystack = A
        b = copy(binarystack)
        b[b .> 0] .= 1
        coms = Images.label_components(b)
        lengths = Images.component_lengths(coms)[2:end]
        indices = Images.component_indices(coms)[2:end]
        boxes = Images.component_boxes(coms)[2:end]
        i = 1
        xr, yr, zr = abs.(boxes[i][1] .- boxes[i][2]) .+ 1
        xy = sqrt(xr^2 + yr^2)
        _xy, _zr, _zp = getextent(boxes[i])
        @test _xy == xy
        @test _zr == zr
        @test _zp == boxes[i][1][3] + zr/2
    end


    @testset "full3dstack" begin
        Random.seed!(42)
        X, Y, Z = 50, 50, 10
        for _ in 1:10
            img3 = Images.N0f8.(rand(X,Y,Z))
            img4 = Images.N0f8.(rand(X,Y,Z))
			@warn "Fix API"
            # compute_contacts_3d(img3, img4, k=3, w=2, deconvolved=true, sigmas=nothing)
            compute_contacts_3d(img3, img4, k=3, w=2, deconvolved=false, sigmas=[1, 1, 1])
        end
    end

    @testset "filter3d" begin
        center1=[50 50 50]
        center2=[50 50 100]
        X, Y, Z = 128, 128, 128
        sigma = 3
        radii = [20, 30]
        _img1, _img2, _img1f, _img2f = makespheres(radii, [center1,center2],  X, Y, Z, sigma)
        _img1f[_img1f .< 0] .= 0
        _img2f[_img2f .< 0] .= 0
        # return rawcontacts, filteredcontacts, gradientcontacts, img_1f, img_2f, sigmap
        # sfp, img_1f, img_2f, zs, ts, _, _, _, _ = process_contact_stack3d(tiffiles, z, _stride)
        rawcontacts, filteredcontacts, gradientcontacts, img1f, img2f, sigmap = compute_contacts_3d(_img1f, _img2f, k=3, w=2, deconvolved=true, sigmas=nothing)
        @test all(rawcontacts .>= 0)
        @test all(rawcontacts[img1f .* img2f .== 0] .== 0)
        @test sum(rawcontacts) >= sum(filteredcontacts)
        @test sum(filteredcontacts) >= sum(gradientcontacts)
        @test sum(_img1f) > sum(img1f)
    end


    @testset "sp3d" begin
        X, Y,Z = 10, 10, 5
        img1 = zeros(X,Y,Z)
        img2 = zeros(X,Y,Z)
        Random.seed!(42)
		for _ in 1:10
	        img3 = rand(X,Y,Z)
	        img4 = rand(X,Y,Z)
	        for stride in 1:3
	            rsp, _zs, _ts = sp3d(img1, img1, stride)
	            @test all(rsp[1+stride:X-stride, 1+stride:Y-stride, 1+stride:Z-stride] .== 1)
	            @test all(rsp[1:1+stride-1, 1:1+stride-1, 1:1+stride-1] .== 0)
	            rrsp, _rzs, _rts = sp3d(img3, img4, stride)
	            @test all(rrsp .<= 1.0)
	            @test all(rrsp .>= -1.0)
	        end
		end
        s0, zns, _ = sp3d(img1, img2, 1, subsampleincrement=0)
        s, zs, _ = sp3d(img1, img2, 1, subsampleincrement=4)
        Statistics.mean(zns[(zns .!= Inf) .& (zns .!= 0)]) < Statistics.mean(zs[(zs .!= Inf) .& (zs .!= 0)])
        @test length(unique(s)) == 2 # 0, 1
        @test maximum(s) == 1
        @test minimum(s) == 0
    end

    @testset "spearsubsample" begin
        X, Y,Z = 10, 10, 5
        img1 = zeros(X,Y,Z)
        img2 = zeros(X,Y,Z)
        Random.seed!(42)
        img3 = rand(X,Y,Z)
        img4 = rand(X,Y,Z)
        expected = [0.011164835, -0.00015337618, 0]
        s, _, _ = spear(img3, img4)
        s2, _, _ = spearsubsample(img3, img4, increment=0)
        @test s==s2
        s2, z1, t1 = spearsubsample(img3, img4, increment=1)
        @test s!=s2
        s2, z0, t0 = spearsubsample(img1, img2, increment=0)
        @test s2 == 1
    end

	@testset "nc" begin
		Random.seed!(42)
        b = rand(100,100) .+ 0.5
		cts = normcontacts(b)
		@test all(cts .>= 0)
		@test minnz(cts) < minnz(b)
	end

	@testset "minr" begin
		# function compute_sample_size_for_min_corr(r, α, β)
		r = 0.5
		α = 0.05
		N = compute_sample_size_for_min_corr(r, α, α)
		@test N > 0
		M = compute_sample_size_for_min_corr(r, α/10, α/10)
		@test M > N
	end

	@testset "iof" begin
		m, s = 1 .+ rand(2)
		io = indexofdispersion(m, s)
		@test !isnan(io)
	end

	@testset "gradm" begin
		x, y, z = rand(3)
		@test gradientmagnitude3d(x, y, z) > 0
	end


    @testset "binarize" begin
        Random.seed!(42)
        b = .5 .+ rand(100,100)
        th = .95
        _th, os, ns = binarize(b, th)

        @test all(b[_th.==1] .== th)
        @test all(os[0 .< b.< th] .== 0.0)
        @test all(os[b.>th] .== 1.0)
        @test all(ns[0 .< b.< th] .== 1.0)
        @test all(ns[b.>=th] .== 0.0)
        @test all(ns .* os .== 0.0)
    end

end
