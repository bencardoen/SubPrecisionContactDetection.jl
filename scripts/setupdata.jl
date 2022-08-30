using Images


p = "/tmp/testdata"
if ! isdir(p)
    mkdir(p)
end
A = zeros(100, 100, 10)
Images.save(joinpath(p, "1.tif"), A)
Images.save(joinpath(p, "2.tif"), A)
