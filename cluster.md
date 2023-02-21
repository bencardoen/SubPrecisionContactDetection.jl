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
A helpful step-by-step guide to allow you to run MCSDetect on 100s of cells at the same time.

## Table of Contents
- [0 Requirements ](#req)
- [1 Login to cluster ](#log)
- [2 Validate dataset ](#val)
- [3 Schedule dataset ](#sched)
- [4 PostProcess results ](#post)
- [5 Export output ](#export)
- [6 Troubleshooting ](#issues)

<a name="req"></a>
## 0. What you will need
- Access to the [cluster](https://ccdb.computecanada.ca/security/login)
- [SSH setup](https://docs.alliancecan.ca/wiki/SSH)
- Data is stored on the cluster in the following format (do NOT use spaces in naming please)
```
- Topdirectory
--- Replicate number (directory named 1, 2, ...)
    --- Treatment (directory named e.g. "GP78+", "HT-1080", ...)
        --- SeriesXYZ (directory that holds the image data for 1 (segmented deconvolved cell), where XYZ is an integer, e.g. "Series123")
            --- channel01.tif (3D 16bit image file containing mitochondria channel, name ending with "1.tif")
            --- channel02.tif (ER channel, ending with "2.tif")
```

Scheduling large (Terabytes) data is time consuming, and you do not want to waste your own time finding out some of your data isn't properly organized.
We will [validate](#val) your dataset to ensure no errors lead to wasted time.

In this guide we will refer to your username on the cluster as `$USER'. 

In Linux command line, $VAR references a variable, e.g. $HOME is your home directory, typically `/home/$USER`.

<a name="log"></a>
## 1. Login
```bash
ssh $USER@cedar.computecanada.ca
```
This should result in something like:
```
[bcardoen@cedar1 ~]$ 
```
(bcardoen is my user name, '~' means you're in your home directory

### 1.1 Viewing variable names
```bash
echo $USER
echo $HOME
```
This will print something like
```
[bcardoen@cedar1 ~]$ echo $USER
bcardoen
[bcardoen@cedar1 ~]$ echo $HOME
/home/bcardoen
[bcardoen@cedar1 ~]$ 
```

We will instruct you to set variables, this is easily done:
```bash
export MYVAR="somevalue"
```
Let's test this, it should show:
```
[bcardoen@cedar1 ~]$ export MYVAR="somevalue"
[bcardoen@cedar1 ~]$ echo $MYVAR
somevalue
```

Great, now let's move on

**OPTIONAL** (but HIGHLY recommended)

Once login is succesful, we will start a `tmux` session to ensure any network interruptions do not break your workflow
```bash
tmux
```

If you want to reconnect to an existing session:
```bash
tmux -t 0  # to reconnect to session 0
```

If you want to view which sessions are active
```bash
tmux list
```

Note that there are multiple login servers, if you can't find your session, ensure you login to the right now (cedar1 vs cedar5).
```bash
ssh $USER@cedar5.computecanada.ca 
```

<a name="val"></a>
## 2. Validate data
### 2.1 Create a new clean working directory
This directory will hold intermediate files needed during processing.

```bash
export EXPERIMENT="/scratch/$USER/myexperiment"
mkdir -p $EXPERIMENT
cd $EXPERIMENT
```

### 2.2 Get DataCurator
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
Alternatively, you can download the latest version
```bash
module load singularity
singularity pull --arch amd64 library://bcvcsert/datacurator/datacurator:latest
```
### 2.3 Acquire computational resource
You'll need your group id, which is of the form `rrg-yourpi` or `def-yourpi`
```bash
export MYGROUP="rrg-mypi" # Replace this with something valid for you
salloc --mem=64GB --account=$MYGROUP --cpus-per-task=16 --time=3:00:00 
```
This will log you in to a compute node with 16 cores, 64GB, for 3 hours.
### 2.4 Copy recipe
DataCurator needs a recipe to verify, this recipe can be found [online](https://github.com/bencardoen/SubPrecisionContactDetection.jl/blob/main/recipe.toml).
```bash
wget https://raw.githubusercontent.com/bencardoen/SubPrecisionContactDetection.jl/main/recipe.toml
```


### 2.5 Update recipe
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
**NOTE** If your channels are 0.tif an 1.tif, rather than 1.tif and 2.tif, please edit the template to reflect this:
```
sed -i "s|[1,2].tif|[0.1].tif|" recipe.toml ## Optional if you need to change channels
```
### 2.6 Configure Slack/Owncloud uploading [Optional]
See [DataCurator documentation](https://github.com/bencardoen/DataCurator.jl/blob/main/docs/src/remote.md)

### 2.7 Validate your data with DataCurator
```bash
module load singularity
export SINGULARITY_BINDPATH="/scratch/bcardoen,$SLURM_TMPDIR"  
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
./datacurator_latest.sif -r recipe.toml
```

This will do the following:
- Build lists for batch processing of all valid files in `in.txt` and `out.txt`
- Report any data that isn't matching the recipe in `errors.txt`
- Compute intensity statistics of all valid data in `channels.csv`
- Compute object statistics of all valid data in `objects.csv`

<a name="sched"></a>
## 3 Schedule the dataset
In your current directory you should now have 2 files
- in.txt
- out.txt

We will ask the cluster to process all cells listed in `in.txt`, with output to be stored in `out.txt`.

### 3.1 Create a scheduling task
Let's get a prepared script to do just that
```bash
wget https://raw.githubusercontent.com/bencardoen/SubPrecisionContactDetection.jl/main/hpcscripts/arraysbatch.sh
```
Make it executable
```
chmod u+x arraysbatch.ch

```
### 3.2 Update with your information (email, group)
You need to change this script to match your account in **3** places:
- EMAIL
- ACCOUNT
- Nr of cells
```bash
export MYEMAIL="you@ubc.ca"
export MYACCOUNT="rrg-mypi"
NCELLS=`wc -l in.txt | awk '{print $1}'`
sed -i "s|CELLS|${NCELLS}|" arraysbatch.sh
sed -i "s|EMAIL|${MYEMAIL}|" arraysbatch.sh
sed -i "s|MYACCOUNT|${MYACCOUNT}|" arraysbatch.sh
sed -i "s|[1,2].tif|[0.1].tif|" arraysbatch.sh ## Optional if you need to change channels
```

### 3.3 Download MCSDetect
Next, we need to make sure the MCS detect singularity image is in place
```bash
singularity pull --arch amd64 mcsdetect.sif library://bcvcsert/subprecisioncontactdetection/mcsdetect_f35_j1.7:j1.8
```
This should download the file **mcsdetect.sif** in your current directory.

### 3.4
Now it's time to submit the job
```bash
sbatch arraysbatch.sh
```
You will get email updates on job progression, and in the current direcotry `.out` files will be saved that contain logs of all the jobs.

Output will be saved in $OUTPUT.

### 3.5 View progress
```bash
squeue -u $USER
```
Will show you the status of your current running jobs.

<a name="post"></a>
## 4 Post

You can now run the postprocessing on the output, see [instructions](https://github.com/bencardoen/SubPrecisionContactDetection.jl/blob/main/scripts/run_cube_sampling_on_dataset.jl).

You can do this with the singularity image (mcsdetect.sif).
Let's create a new directory where to store the postprocessing:
```bash
export POSTDATA="/scratch/$USER/post"
mkdir -p $POSTDATA
```
We'll need the data output directory you used earlier:
```
export OUTPUT="/project/myresearchgroup/myoutput"
```
Navigate to the directory you created earlier
```bash
cd $MYEXPERIMENT
```
**make sure the mcsdetect.sif is in this directory**

Download a script for you to do the postprocessing
```bash
wget https://raw.githubusercontent.com/bencardoen/SubPrecisionContactDetection.jl/main/hpcscripts/postprocess.sh
chmod u+x postprocess.sh
```
Update the script with your variables:
**Please make sure the variables are correctly set as done above**
```
sed -i "s|EMAIL|${MYEMAIL}|" postprocess.sh        # Your email
sed -i "s|MYACCOUNT|${MYACCOUNT}|" postprocess.sh  # Your cluster account
sed -i "s|POSTOUTPUT|${POSTDATA}|" postprocess.sh  # The location where you want the postprocessed results
sed -i "s|INPUT|${OUTPUT}|" postprocess.sh         # The computed contacts location, previously saved in $OUTPUT
```
Then submit
```bash
sbatch postprocess.sh
```
**If you want to change the parameters, please edit the script (using nano, vim, ..)**

You'll get an email when your results are done.

<a name="export"></a>
## 5 Export your results
I recommend using [Globus](https://globus.computecanada.ca/) to export results efficiently.

- Log in to [Globus](https://globus.computecanada.ca/) using your cluster ID
- Select `$POSTDATA` as the source for the transfer
- Select a directory on your target computer as the target
- Execute
You'll get an email if the transfer completed. 

### 5.1 [Optional] Compress the results
See [wiki](https://docs.alliancecan.ca/wiki/Archiving_and_compressing_files)

### 5.2 [Optional] Selectively extract results
You can use [DataCurator](https://github.com/bencardoen/DataCurator.jl) to write a recipe that, for example, only collects CSV files of the contacts, and sends them to OwnCloud. Or compute statistics of the contacts, mitochondria, etc, etc. 
Examples are in DC documentation, and execution is the same as in step 3.

<a name="issues"></a>
## 6 Troubleshooting

If you run into problems or something is not clear, please create a [new issue](https://github.com/bencardoen/SubPrecisionContactDetection.jl/issues/new/choose)
