# -*- coding: utf-8 -*-
"""
Created on Thu May 20 11:41:50 2021

@author: kbefus
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

wdir = os.path.join(r'O:\seepage_face_UZF','models','gen_seep')
model_name_fmt = 'rel{0:3.2f}_vs{1:.0f}_aqb{2:.0f}_RK{3:3.2f}'

vdepth = 50.
vwidth = 500.
rel_lev = 0.4
aq_thick = 100.
rdK = 0.03
log10rdK = np.log10(rdK)
L=5e3
nx=1000


model_name = model_name_fmt.format(rel_lev,vdepth/vwidth,L,log10rdK)
mdir = os.path.join(wdir,model_name)

if not os.path.isdir(mdir):
    os.makedirs(mdir)

in_dict={'slope':0.1,'depth':0.12*500,'width':500,
         'rel_lev':0,'aq_thick':aq_thick,
         'r_over_K':rdK,'wdir':mdir,'ztype':'spread',
         'Lx':L,'ncol':nx,'use_drn':False,'kmult':1}
#%%
icount=0
in_dict['modelname'] = '{}_{}'.format(model_name,icount)
m_out,elev,cell_types = seep_model(**in_dict)

cbc_fname = os.path.join(mdir,'{}.cbc'.format(in_dict['modelname']))
cbb = bf.CellBudgetFile(cbc_fname)
rch_in = cbb.get_data(text='UZF RECHARGE')[0].squeeze()
dr = cbb.get_data(text='SURFACE LEAKAGE')[0].squeeze()
cbb.close()
cbb = None

plt.close('all')
plt.plot(rch_in.sum(axis=0),'b.',label='recharge flux')
plt.plot(dr[0],'rx',label='drain flux')
plt.legend()
#%%
output=seep_face_properties(modelname=in_dict['modelname'] ,wdir=mdir,dem=elev,cell_types=cell_types,totim=1,dx=L/nx)