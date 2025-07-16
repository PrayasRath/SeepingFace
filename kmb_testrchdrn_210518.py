# -*- coding: utf-8 -*-
"""
Created on Tue May 18 11:21:57 2021

@author: Kbefus
"""

import os,sys
import matplotlib.pyplot as plt
import numpy as np

lib_path = r'E:\seepage_face_UZF'
sys.path.insert(0,lib_path)

from seepfuncs_rath_kmb_210517 import seep_model,seep_test,plot_model

#%%

wdir =  os.path.join(r'O:\seepage_face_UZF','models','gen_seep')
model_name_fmt = 'rel{0:3.2f}_vs{1:.0f}_aqb{2:.0f}_RK{3:3.2f}'

vdepth = 50.
vwidth = 500.
rel_lev = 0.6
aq_thick = 100.
rdK = 1e-1
log10rdK = np.log10(rdK)
L=1e4
nx=1000


model_name = model_name_fmt.format(rel_lev,vdepth/vwidth,aq_thick,log10rdK)
mdir = os.path.join(wdir,model_name)

if not os.path.isdir(mdir):
    os.makedirs(mdir)

in_dict={'slope':vdepth/vwidth,'depth':vdepth,'width':vwidth,
         'rel_lev':rel_lev,'aq_thick':aq_thick,
         'r_over_K':rdK,'wdir':mdir,'ztype':'spread',
         'Lx':L,'ncol':nx,'use_drn':True}
#%%
icount = 0
max_iter = 20
ndrain_shutoff = 5
exit_bool = False
drn_bool_last = None
force_run = True
all_drn_bools = []
all_drn_used = []
while not exit_bool:
    print('Model iteration = {}'.format(icount))
    in_dict['modelname'] = '{}_{}'.format(model_name,icount)
    m_out,elev,cell_types = seep_model(**in_dict)
    # all_ct.append(cell_types)
    drn_bool = seep_test(in_dict['modelname'],mdir,cell_types)
    all_drn_bools.append(drn_bool)
    if drn_bool_last is not None:
        active_drn = drn_bool.nonzero()[0]
        diff_list1 = list(set(active_drn)-set(drn_bool_last.nonzero()[0]))
        diff_list2 = list(set(drn_bool_last.nonzero()[0])-set(active_drn))
        if len(diff_list2) == 0 and len(diff_list1)==0 and not force_run: # both datasets have exact same drains
            exit_bool = True
        else:
            icount += 1
            # remove ndrain_shutoff rch cell(s) per loop
            # in_dict['drn_bool'] = drn_bool.copy() # use in next model run
            in_dict['drn_bool'] = drn_bool_last.copy()
            in_dict['drn_bool'][drn_bool_last.nonzero()[0][-ndrain_shutoff]:] = False  #turns off recharge to cells where True
            drn_bool_last = in_dict['drn_bool'].copy()
    else:
        icount += 1
        # in_dict['drn_bool'] = drn_bool.copy() # use in next model run
        in_dict['drn_bool'] = np.ones_like(drn_bool,dtype=bool) # all True, no recharge cells
        in_dict['drn_bool'][in_dict['drn_bool'].nonzero()[0][-ndrain_shutoff]:] = False
        drn_bool_last = in_dict['drn_bool'].copy()
    
    all_drn_used.append(in_dict['drn_bool'])
    
    if icount>max_iter:
        exit_bool = True
        print("max iterations reached")
ad = np.array(all_drn_bools)
adu = np.array(all_drn_used)
#%%
plt.close('all')
fig1,ax1 = plt.subplots(3,1)
ax1 = np.array(ax1).ravel()
c1=ax1[0].matshow(ad.astype(int)) # model output of drain locations (1=Drain, 0=No drain)
plt.colorbar(c1,ax=ax1[0],location='bottom')
ax1[0].set_title('Model predictions of active drain locations, 1=Drain, 0=No drain')
ax1[0].set_ylabel('model iteration')


c2= ax1[1].matshow(adu.astype(int)) # model prescribed locations of drains (1=Drain/no rch, 0=Rch/no drain)
ax1[1].set_title('Prescribed recharge locations, 1=No Rch, 0 = Rch')
plt.colorbar(c2,ax=ax1[1],location='bottom')
ax1[1].set_ylabel('model iteration')

ax1[2].set_title('0=Both rch or both drn; 1=No drn and no rch; -1: drn incorrectly w/ rch')
c3=ax1[2].matshow(adu.astype(int)-ad.astype(int)) # where the two are different
plt.colorbar(c3,ax=ax1[2],location='bottom')
ax1[2].set_xlabel('model column')
ax1[2].set_ylabel('model iteration')
#%%
plot_iter = 0
fig2,ax2 = plt.subplots()
plt_modelname = in_dict['modelname'].replace('_{0:.0f}'.format(icount-1),'_{0:.0f}'.format(plot_iter))
plot_model(plt_modelname,mdir,ax=ax2)
ax2.set_title(plt_modelname)
# xtemp = np.arange(len(elev)) * L/nx + (L/nx)/2
# ax2.plot(xtemp[adu[plot_iter]],1+elev[adu[plot_iter]],'ko')
# ax2.plot(xtemp[~adu[plot_iter]],1+elev[~adu[plot_iter]],'bx')