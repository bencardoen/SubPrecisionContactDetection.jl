## This file documents how you can run interactive Julia commands
## To execute a line, place the cursor at the end of the line or highlight it, and click CTRL+ENTER. 
## A ↺ will appear, followed by a ✓

## Start the session
using Pkg; Pkg.activate(".");

## Load the packages we need
using SubPrecisionContactDetection, Images, DataFrames

## Let's say you want to test a few alpha values
# Get the default parameters
parameters = get_defaults();
​
# Let's check what the values are:
for key in keys(parameters)
    @info "Key $key is set to $(parameters[key])"
end
​
inpath = "C:\\yourdata"
outpath = "C:\\youroutput"
# Define your values
alphas = [0.05, 0.01, 0.1]
​
# Run in a simple for loop
for alpha in alphas
    # Get the default values (30)
    params = get_defaults()
    # Set the input
    params["inpath"]=inpath
    # Let's create a subfolder for each alpha
    op = mkpath(joinpath(outpath, "$alpha"))
    # Set the output to this path
    params["outpath"]=outpath
    # 2D, so let's tell to code
    params["dimension"]=2
    # And it's deconvolved, so let's set that to true
    params["deconvolved"]=true
    two_channel_contacts(params)
end

### Loading and saving images
img = Images.load("my.tif");

a = zeros(10, 10, 10);
# Make sure it's [0-1]
clamp01!(a);
Images.save("1.tif", Images.N0f16.(a));

### Computing object descriptors (features)
dataframe = describe_objects(img);

### Save to CSV
CSV.write("test.csv", dataframe)