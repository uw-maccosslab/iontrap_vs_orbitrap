#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 16:56:55 2020

@author: Lilian
"""

from pyteomics import mzml
import numpy as np
import pandas as pd

# Specify file name like so:
filename="quant_comp_lastscanagc/20201024_LRH_plasma_OT_0-7_5.mzML"



mzml_file=mzml.read(filename, use_index=True)



df=np.array([[mzml_file[0]['id'].split("scan=",1)[1],
              mzml_file[0]['ms level'],
              mzml_file[0]['total ion current'],
              mzml_file[0]['scanList']['scan'][0]['ion injection time'],
              mzml_file[0]['scanList']['scan'][0]['scan start time']]])

for x in range(1, len(mzml_file)):
    if mzml_file[x]['ms level']==2:
        df=np.append(df, [[mzml_file[x]['id'].split("scan=",1)[1],
                           mzml_file[x]['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z'],
                           mzml_file[x]['total ion current'],
                           mzml_file[x]['scanList']['scan'][0]['ion injection time'],
                           mzml_file[x]['scanList']['scan'][0]['scan start time']]], axis = 0)

df=pd.DataFrame(df, columns = ['scan', 'precursormz', 'totalIonCurrent',
                               'injectTime', 'scanTime'])
df=df.drop(df.index[0])

# Change output name and save file
df.to_csv("quant_plasma_OT_0-7_5.csv")
