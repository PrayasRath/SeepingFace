# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:38:11 2021

@author: kbefus
"""
import os, sys
import matplotlib.pyplot as plt
import numpy as np

lib_path = r'E:\seepage_face_UZF'
sys.path.insert(0,lib_path)

from seepfuncs_rath_kmb_210517 import seep_active_regional,seep_active
#%%
wdir = os.path.join(r'E:\seepage_face_UZF','models','gen_seep_func')
in_dict = {'slope':0.1,
           'depth':50,
           'width':500,
           'aq_thickness':100,
           'nx':1000, # use same number of cells per model
           'Lx':5e3,
           'ztype':'spread', # zbot type = ['top','bot','spread']
           'wdir':wdir,
           'rel_level':[0.5,0.75],# not all of them
           'RK_array':[1E-5,1E-4,1E-3,1E-2,1E-1,1][::-1],
           }

# Would be better to rely mainly on a dataframe for outputs
all_RK_list_regio,seep_active_all_RL,data_df = seep_active_regional(**in_dict)

# To select only the onset of seepage models
onset_df = data_df.loc[data_df['onset_flag'],:].copy()
