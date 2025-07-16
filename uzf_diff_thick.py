# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 15:05:14 2021

@author: prath1
"""
import os,sys
import matplotlib.pyplot as plt
import numpy as np
from flopy.utils import binaryfile as bf
import flopy.utils.postprocessing as pp
import matplotlib.pyplot as plt
from flopy.plot import PlotCrossSection as pcs

lib_path = r'E:\seepage_face_UZF'
sys.path.insert(0,lib_path)

from seepfuncs_rath_kmb_210517 import seep_model,seep_test,plot_model,bf,seep_face_properties

#%%

wdir = os.path.join(r'E:\seepage_face_UZF2','SLF_SFF','diff_thick_hi_slope')
model_name_fmt = 'rel{0:3.2f}_vs_thickness{1:3.2f}'

slope=0.1 
vdepth = 50.
width=500
rel_lev = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
thickness_list=np.array([0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15])
rdK =0.03
L=5e3
nx=1000
SLF_list=[]
SFF_list=[]
SLF_regional_list=[]
SFF_regional_list=[]
for i in range (len(thickness_list)):
    SLF=[]
    SFF=[]
    SLF_regional=[]
    SFF_regional=[]
    for j in range (len(rel_lev)):
        
        model_name = model_name_fmt.format(rel_lev[j],thickness_list[i])
        mdir = os.path.join(wdir,model_name)

        if not os.path.isdir(mdir):
            os.makedirs(mdir)

        in_dict={'slope':0.1,'depth':50,'width':width,
         'rel_lev':rel_lev[j],'aq_thick':thickness_list[i]*L,
         'r_over_K':rdK,'wdir':mdir,'ztype':'spread',
         'Lx':L,'ncol':nx,'use_drn':False,'kmult':1}

        in_dict['modelname'] = model_name
        m_out,elev,cell_types = seep_model(**in_dict)
        all_out=seep_face_properties(modelname=model_name,wdir=mdir,cell_types=cell_types,dem=elev,dx=L/nx)
        dr_seep_total,water_table,sff,slf,slf_2d,dr_seep_out_total,ch_ratio,ch_out,q_ratio,sff_regional,slf_regional=all_out
        SLF.append(slf)
        SFF.append(sff)
        SLF_regional.append(slf_regional)
        SFF_regional.append(sff_regional)
        print(model_name)
    SLF_list.append(SLF)
    SFF_list.append(SFF)
    SLF_regional_list.append(SLF_regional)
    SFF_regional_list.append(SFF_regional)
os.chdir(wdir)
np.savetxt("uzf_diff_thick_SLF_hi_s.csv",SLF_list , delimiter=",", fmt='%s')
np.savetxt("uzf_diff_thick_SFF_hi_s.csv",SFF_list , delimiter=",", fmt='%s')
np.savetxt("uzf_diff_thick_SLF_regional_hi_s.csv",SLF_regional_list , delimiter=",", fmt='%s')
np.savetxt("uzf_diff_thick_SFF_regional_hi_s.csv",SFF_regional_list , delimiter=",", fmt='%s')
