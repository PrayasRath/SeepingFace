# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 11:58:50 2021

@author: prath1
"""
import os,sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
lib_path = r'E:\seepage_face_UZF2'
sys.path.insert(0,lib_path)

from seepfuncs_rath_kmb_210517 import seep_active_regional,seep_active
#%%
wdir = os.path.join(r'E:\seepage_face_UZF2','onset_of_gwd','slope_ratio_effect7')
seep_active_all_list=[]
regional_onset_df_diff_slope=[]
regional_slope_list=np.array([10**(3),10**(2.5),10**(2),10**(1.5),10**(1),10**(0.5),
                              10**(0.),10**(-0.5),10**(-1),10**(-1.5),10**(-2)])
for i in range (len(regional_slope_list)):
     newpath=os.path.join(wdir,"slope_ratio"+str(regional_slope_list[i]))
     if not os.path.exists(newpath):
        os.makedirs(newpath)
     else:
        os.chdir(newpath) 
     in_dict = {'slope':(5/1000)/regional_slope_list[i],
               'depth':5,
               'width':1000,
               'aq_thickness':100,
               'nx':1000, # use same number of cells per model
               'Lx':5e3,
               'ztype':'spread', # zbot type = ['top','bot','spread']
               'wdir':newpath,
               'rel_level':[0.,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99],# not all of them
               'RK_array':[1E-5,1E-4,1E-3,1E-2,1E-1,1][::-1],
               }

# Would be better to rely mainly on a dataframe for outputs
     all_RK_list,seep_active_all_RL,data_df = seep_active_regional(**in_dict)
     seep_active_all_list.append(seep_active_all_RL)
# To select only the onset of seepage models
     regional_onset_df_diff_slope.append( data_df.loc[data_df['onset_flag'],:].copy())
os.chdir(newpath)
np.savetxt("uzf_regional_onset_ratio_slope7.csv",seep_active_all_list , delimiter=",", fmt='%s')
print( pd.concat(regional_onset_df_diff_slope).to_csv('valley_regional_diff_reg_slope7.csv'))