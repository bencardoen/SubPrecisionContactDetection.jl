<!--
 Copyright (C) 2023 Ben Cardoen
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.
 
 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
-->

# MCSDetect on large datasets on clusters
A helpful guide to allow you to run MCSDetect on 100s of cells at the same time.

In this guide we will refer to your username on the cluster as `$USER'. 

In Linux command line, $VAR references a variable, e.g. $HOME is your home directory, typicall `/home/$USER`.


## Login
```bash
ssh $USER@cedar.computecanada.ca
```

## Validate data
### Create a new clean working directory
This directory will hold intermediate files needed during processing.
```bash
export EXPERIMENT="/scratch/$USER/myexperiment"
mkdir -p $EXPERIMENT
cd $EXPERIMENT
```

### Get DataCurator
For the next step we'll need to download DataCurator to validate your dataset layout.
You can obtain it [here](https://github.com/bencardoen/DataCurator.jl), but it will be present on Cedar.
```bash
module load singularity
# Let's copy the executable image to our current location  "."
cp /project/rrg-hamarneh/singularity_images/datacurator_latest.sif . 
# Make sure it's here
ls -lahst *.sif
```
Should show something like this
```bash
2.3G -rwxr-x--- 1 bcardoen rrg-hamarneh 2.5G Feb 19 07:27 datacurator_latest.sif
```
Make sure it's executable
```bash
chmod u+x datacurator_latest.sif
```
### Acquire computational resource
You'll need your group id, which is of the form `rrg-yourpi` or `def-yourpi`
```bash
export MYGROUP="rrg-mypi" # Replace this with something valid for you
salloc --mem=64GB --account=$MYGROUP --cpus-per-task=16 --time=3:00:00 
```
This will log you in to a compute node with 16 cores, 64GB, for 3 hours.
### Copy recipe
DataCurator needs a recipe to verify, this recipe can be found online, for your convenience this is what it should look like:
```toml
[global]
# act_on_success=false
inputdirectory = "INPUT"
hierarchical=true
traversal="topdown"
parallel=true
regex=true
common_actions = {do_on_invalid=[["all", "show_warning", ["log_to_file", "errors.txt"]]]}
common_conditions = {is_3d_channel=[["all", "is_tif_file", "is_3d_img", "filename_ends_with_integer"]], dir_only=[["all", "isdir", ["not", "is_hidden"]]]}
file_lists=[{name="in", aggregator=[["filepath","sort","unique","shared_list_to_file"]]},
              {name="out", aggregator=[[["change_path", "OUTPUT"],"filepath","sort","unique","shared_list_to_file"]]},
              {name="objects", aggregator=[["describe_objects","concat_to_table"]]},
              {name="channels", aggregator=[["describe_image","concat_to_table"]]},]
[any]
conditions=["never"]
actions=["do_on_invalid"]

[level_1]
conditions=["dir_only"]
actions=["do_on_invalid"]

[level_2]
all=true
conditions=["dir_only", "ends_with_integer", ["startswith", "Replicate"]]
actions=["do_on_invalid"]

[level_3]
conditions=["dir_only"]
actions=["do_on_invalid"]

[level_4]
all=true
conditions=["dir_only", ["startswith", "Serie"], "ends_with_integer"]
actions=["do_on_invalid"]

[level_5]
conditions=["is_3d_channel"]
actions=["do_on_invalid"]
counter_actions=[["->", ["in","out","objects","channels"]]]
```

### Update recipe
We need to update this template with the data locations:
Let us assume the data you want to process is located in
```bash
export DATA="/project/myresearchgroup/mydata"
export OUTPUT="/project/myresearchgroup/myoutput"
```
You can either do this with an editor, or with these commands
```bash
sed -i "s|INPUT|${DATA}|" recipe.toml # Set the correct data directory
sed -i "s|OUTPUT|${OUTPUT}|" recipe.toml
```
If you now check the recipe these 2 lines should be changed to your data
```bash
inputdirectory = "INPUT"
# will be
inputdirectory = "/project/myresearchgroup/mydata"
# and
              {name="out", aggregator=[[["change_path", "OUTPUT"],"filepath","sort","unique","shared_list_to_file"]]},
# will be
              {name="out", aggregator=[[["change_path", "/project/myresearchgroup/myoutput"],"filepath","sort","unique","shared_list_to_file"]]},
```

### Run DataCurator
```bash
module load singularity
export SINGULARITY_BINDPATH="/scratch/bcardoen,$SLURM_TMPDIR"  
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
./datacurator_latest.sif -r recipe.toml
```

This will do the following:
- Build lists for batch processing of all valid files in `in.txt` and `out.txt`
- Compute intensity statistics of all valid data in `channels.csv`
- Compute object statistics of all valid data in `objects.csv`