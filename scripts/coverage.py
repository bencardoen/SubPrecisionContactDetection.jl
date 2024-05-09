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
# Copyright 2021-2024, Ben Cardoen
import pandas as pd
import numpy as np
import argparse
import logging
import os
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
lgr = logging.getLogger('global')
lgr.setLevel(logging.INFO)
lgr = None

def postprocess_sampled(_df):
#     CUBEDF =  pd.read_csv(os.path.join(path, "all.csv"))
    CUBEDF=_df.copy()
    CUBEDF['ratio_cf_to_mf'] = CUBEDF['ctsurface']/CUBEDF['mtsurface'] * 100
    CUBEDF['mean_mito'] = CUBEDF['mitosum']/CUBEDF['mitvol']
    CUBEDF['mean_spear']=CUBEDF['contactsum'] / CUBEDF['contactvol']
    CUBEDF['mean_spear'].fillna(0, inplace=True)
    return CUBEDF.copy()

nq = lambda x : np.quantile(x, .75)
def aggregate(df):
    ### Data is organized by >Replicate>Celltype>Serienr
    q = df.groupby(
        ['serienr', 'celltype', 'replicate', 'experiment']
    ).agg(
        {
            'mitvol' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt],
            'contactsum' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt],
            'ncontacts' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt, nq],
            'ctsurface' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt, nq],
            'ratio_cf_to_mf' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt, nq],
            'mitosum' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt, nq],
            'mtsurface' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt, nq],
            'contactvol' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt, nq],
            'mean_mito' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt, nq],
            'mean_spear' : ['mean', 'std', 'count', 'sum', 'skew', pd.Series.kurt, nq],
        }
        ).reset_index()
    q.columns = [' '.join(col).strip() for col in q.columns.values]
    return q





# From https://github.com/bencardoen/ERGO.py/blob/main/src/gconf.py
def initlogger(configuration):
    global lgr
    if lgr is None:
        lgr = logging.getLogger('global')
    if 'logdir' in configuration:
        fh = logging.FileHandler(os.path.join(configuration['logdir'], 'svrg.log'))
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter("[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s")
        fh.setFormatter(formatter)
        lgr.addHandler(fh)
    lgr.setLevel(logging.INFO)
    return lgr

# From https://github.com/bencardoen/ERGO.py/blob/main/src/gconf.py
def getlogger():
    global lgr
    if lgr is None:
        return initlogger({})
    return lgr


import os
def run(args):
    getlogger().info("Loading input CSV file ...")
    _ALLDF = pd.read_csv(os.path.join(args.inputdirectory, "all.csv"))
    getlogger().info("... Done")
    ALLDF=_ALLDF.copy()
    ALLDF['experiment'] = args.inputdirectory
    ALLDF = postprocess_sampled(ALLDF.copy())
    getlogger().info("Describing data")
    aggregated = aggregate(ALLDF)
    EXPS = np.unique(aggregated['experiment'])
    for exp in EXPS:
        _DF = aggregated[aggregated['experiment'] == exp]
        rs = np.unique(_DF['replicate'])
        print("Have {} replicates for {}".format(rs, exp))
        _CD = _DF[_DF['replicate']<99].copy()
        for ct in np.unique(_CD['celltype']):
            ccount = len(_CD[_CD['celltype'] == ct])
            print("Experiment {} CT {} has {} cells".format(exp, ct, ccount))
    

    FDF = ALLDF.copy()
    _ADF = FDF[FDF['ctsurface'] > 0]
    _ADF = _ADF[_ADF['mitvol'] > 5]
    SDF=aggregate(_ADF.copy())
    SDF['coverage'] = SDF['contactvol sum']/SDF['mitvol sum'] * 100
    cols = ['serienr', 'celltype', 'replicate', 'experiment', 'mitvol mean',
        'mitvol std', 'mitvol sum','ncontacts mean', 'ncontacts sum','mean coverage % per sliding window','Coverage % mito by contacts, mean per cell','mitvol sum','contactvol sum']
    df=SDF.copy()
    df = df.rename(columns={'ratio_cf_to_mf mean': 'mean coverage % per sliding window'})
    df = df.rename(columns={'coverage': 'Coverage % mito by contacts, mean per cell'})
    new_df = df.loc[:, cols]
    getlogger().info("Saving to{}".format(os.path.join(args.outputdirectory, "coverage_aggregated.csv")))
    new_df.to_csv(os.path.join(args.outputdirectory, "coverage_aggregated.csv"))
    lgr.info("Done !!")

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script processes aggregate coverage data from MCS detect. It needs a file named all.csv, produced by the run_cube_sampling_on_dataset.jl script. It is used to report per cell coverage and contacts per sliding window')
    # Add arguments
    parser.add_argument('--inputdirectory', required=True, help='where to find all.csv, the output from the sampling code')
    parser.add_argument('--outputdirectory', required=True, help='path to outputdirectory to save postprocess csvs')
    # parser.add_argument('--lnsize', type=float, default=9, help='Minimum size of adjacent mitochondria (natural log, default 9)')
    # parser.add_argument('--mitoint', type=float, default=0.2, help='Minimum intensity (mean) of adjacent mitochondria (default 0.2)')
    parser.add_argument('--alpha', type=float, default=0.05, help='Alpha value to load (0.05 is default)')
    args = parser.parse_args()
    lgr=getlogger()
    for arg in vars(args):
        lgr.info("{} --> {}".format(arg, getattr(args, arg)))

    if not os.path.exists(args.inputdirectory) or not os.path.exists(args.outputdirectory):
        lgr.error("Input path or output path does not exist")
        exit(-1)

    run(args)
