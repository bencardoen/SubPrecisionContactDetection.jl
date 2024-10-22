    ## Let's say you want to test a few alpha values
    # Get the default parameters
    parameters = get_defaults();
    ​
    # Let's check what the values are:
    for key in keys(parameters)
        @info "Key $key is set to $(parameters[key])"
    end
    ​
    inpath = "C:\yourdata"
    outpath = "C:\youroutput"
    # Define your values
    alphas = [0.05, 0.01, 0.1]
    ​
    # Run in a simple for loop
    for α in alphas
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