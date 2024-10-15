using Pkg; Pkg.activate(".")
using Statistics, Distributions, Images, Glob, CSV, DataFrames
using SubPrecisionContactDetection
using ImageFiltering

@info "Creating 3 spheres"

i = zeros(40, 40, 40)
i[2:4, 2:4, 2:4] .= 1
i[20:28, 20:28, 20:22] .= 1
i[35:37, 35:37, 35:37] .= 1
sigma = 2
img = Images.imfilter(i, Kernel.gaussian((sigma,sigma,sigma)))
img[img .>= 1] .= 1
img[img .< 0] .= 0
if isdir("qrx")
    rm("qrx"; recursive=true)
end
mkdir("qrx")
Images.save(joinpath("qrx", "1.tif"), Images.N0f16.(img)