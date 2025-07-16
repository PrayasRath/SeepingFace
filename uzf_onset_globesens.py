# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 22:10:42 2021

@author: prath1
"""
import os,sys
import numpy as np
# from skimage import io
import flopy
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
import pandas as pd
# import flopy.utils.geometry as geo
import flopy.utils.postprocessing as pp
import math
wdir = os.path.join(r'E:\seepage_face_UZF2','Global_sens')
os.chdir(wdir)
from SALib.sample import saltelli
from SALib.analyze import sobol
lib_path = r'E:\seepage_face_UZF'
sys.path.insert(0,lib_path)
from seepfuncs_rath_kmb_210517 import seep_model,seep_test,plot_model,bf,seep_face_properties,seep_active,seep_active_regional
wdir = os.path.join(r'E:\seepage_face_UZF2','Global_sens_onset')
os.chdir(wdir)
#    "# Define the model inputs\n",
problem = {
        'num_vars': 5,
        'names': ['slope','d/w','w/l','rel_lev','B/L'],
        'bounds': [[-4, -0.5], #log(slope)
                   [5, 150], #depth
                   [100, 2500], #width",
                   [0, 0.99], #rl",
                   [10,750], #aq_thick",
                   
                   ]}  
paramValues=saltelli.sample(problem,1000)
paramValues[55]=np.array([-6.07666016e-01,  5.21240234e+01,  1.20976562e+03,  3.18559570e-01,
        4.43261719e+02])
paramValues[56]=np.array([-6.07666016e-01,  5.81713867e+01,  2.07460938e+03,  3.18559570e-01,
        4.43261719e+02])
paramValues[57]=np.array([-6.07666016e-01,  5.81713867e+01,  1.20976562e+03,  6.05698242e-01,
        4.43261719e+02])
paramValues[59]=np.array([-6.07666016e-01,  5.81713867e+01,  1.20976562e+03,  3.18559570e-01,
        4.43261719e+02])
paramValues[638]=np.array([-6.58935547e-01,  9.21557617e+01,  1.58710938e+03,  2.18979492e-01,
        1.75878906e+02])

paramValues[639]=np.array([-6.58935547e-01,  1.11696777e+02,  1.32226562e+03,  2.18979492e-01,
        1.75878906e+02])

paramValues[640]=np.array([-6.58935547e-01,  1.11696777e+02,  1.58710938e+03,  3.34028320e-01,
        1.75878906e+02])

paramValues[1390]=np.array([-5.25634766e-01,  1.62573242e+01,  1.49101562e+03,  7.74887695e-01,
        2.33691406e+02])
paramValues[2953]=np.array([-6.76025391e-01,  4.99584961e+01,  1.97148438e+03,  1.22299805e-01,
        2.04785156e+02])
paramValues[2962]=np.array([-6.76025391e-01,  3.15502930e+01,  2.00664062e+03,  6.93676758e-01,
        2.04785156e+02])
paramValues[5634]=np.array([-5.42724609e-01,  7.71459961e+01,  1.27070312e+03,  2.74086914e-01,
        2.26464844e+02])
paramValues[7035]=np.array([-6.28173828e-01,  5.37817383e+01,  1.23085938e+03,  3.99770508e-01,
        5.01074219e+02])
paramValues[7036]=np.array([-6.28173828e-01,  5.37817383e+01,  2.24570312e+03,  3.58198242e-01,
        5.01074219e+02])
paramValues[7805] =np.array([-6.69189453e-01,  6.22778320e+01,  1.28007812e+03,  4.96450195e-01,
        1.18066406e+02])

paramValues[7806]=np.array([-6.69189453e-01,  4.89672852e+01,  1.16523438e+03,  9.65346680e-01,
        1.18066406e+02])

paramValues[7807]=np.array([-2.38842773e+00,  6.22778320e+01,  1.16523438e+03,  9.65346680e-01,
        1.18066406e+02])

paramValues[7808]=np.array([-2.38842773e+00,  4.89672852e+01,  1.28007812e+03,  9.65346680e-01,
        1.18066406e+02])

paramValues[7809]=np.array([-2.38842773e+00,  4.89672852e+01,  1.16523438e+03,  4.96450195e-01,
        1.18066406e+02])
paramValues[8189]=np.array([-5.05126953e-01,  2.37622070e+01,  1.99257812e+03,  4.78564453e-02,
        4.64941406e+02])

paramValues[8190]=np.array([-5.05126953e-01,  7.84204102e+01,  2.17773438e+03,  2.17529297e-02,
        4.64941406e+02])
paramValues[9285]=np.array([-6.72607422e-01,  2.14965820e+01,  7.48046875e+02,  6.60805664e-01,
        3.63769531e+02])
paramValues[9287]=np.array([-6.72607422e-01,  2.14965820e+01,  7.48046875e+02,  9.71147461e-01,
        3.63769531e+02])
paramValues[10558]=np.array([-5.35888672e-01,  4.75512695e+01,  2.34179688e+03,  5.30288086e-01,
        5.15527344e+02])
paramValues[11615]=np.array([-7.68310547e-01,  2.88598633e+01,  1.41367188e+03,  6.11499023e-01,
        9.77050781e+01])
paramValues[11671]=np.array([-5.49560547e-01,  5.74633789e+01,  6.63671875e+02,  1.78374023e-01,
        7.39550781e+02])

paramValues[11672]=np.array([-5.49560547e-01,  9.22973633e+01,  2.22226562e+03,  1.78374023e-01,
        7.39550781e+02])

paramValues[11673] =np.array([-5.49560547e-01,  9.22973633e+01,  6.63671875e+02,  7.34282227e-01,
        7.39550781e+02])

paramValues[11674]=np.array([-5.49560547e-01,  9.22973633e+01,  6.63671875e+02,  1.78374023e-01,
        4.18662109e+02])

paramValues[11675]=np.array([-5.49560547e-01,  9.22973633e+01,  6.63671875e+02,  1.78374023e-01,
        7.39550781e+02])
# paramValues[12024]=np.array([-6.07666016e-01,  3.70727539e+01,  1.80976562e+03,  5.64125977e-01,
#         3.70996094e+02])

# paramValues[12025]=np.array([-7.13623047e-01,  3.70727539e+01,  1.80976562e+03,  5.64125977e-01,
#         3.70996094e+02])

# paramValues[12026]=np.array([-6.07666016e-01,  7.19067383e+01,  1.80976562e+03,  5.64125977e-01,
#         3.70996094e+02])

# paramValues[12027]=np.array([-6.07666016e-01,  3.70727539e+01,  1.75117188e+03,  5.64125977e-01,
#         3.70996094e+02])

# paramValues[12028]=np.array([-6.07666016e-01,  3.70727539e+01,  1.80976562e+03,  9.36342773e-01,
#         3.70996094e+02])

# paramValues[12332]=np.array([-5.50415039e-01,  3.05944824e+01,  1.26191406e+03,  6.00139160e-01,
#         4.81201172e+02])

# paramValues[12333]=np.array([-5.50415039e-01,  3.05944824e+01,  5.14257812e+02,  4.91857910e-01,
#         4.81201172e+02])
# paramValues[12334]=np.array([ -0.55041504,  30.59448242, 514.2578125 ,   0.60013916,
#         161.48925781])

# paramValues[12335]=np.array([ -0.55041504,  30.59448242, 514.2578125 ,   0.60013916,
#         148.12011719])
# paramValues[12411]=np.array([-7.16186523e-01,  4.14270020e+01,  9.64257812e+02,  5.53732910e-01,
#         1.52392578e+02])

# paramValues[12412]=np.array([ -0.71618652,  41.42700195, 211.9140625 ,   0.90951416,
#         152.3925781])
# paramValues[13949]=np.array([-6.75170898e-01,  1.70715332e+01,  7.08789062e+02,  7.35490723e-01,
#         4.52294922e+02])

# paramValues[13950]=np.array([-6.75170898e-01,  5.10620117e+00,  1.02988281e+03,  9.28850098e-01,
#         4.52294922e+02])

# paramValues[13951]=np.array([-2.23205566e+00,  1.70715332e+01,  1.02988281e+03,  9.28850098e-01,
#         45.2294922e+01])

# paramValues[13952]=np.array([ -2.23205566,   5.10620117, 708.7890625 ,   0.9288501 ,
#         452.2949219])

# paramValues[13953]=np.array([-2.23205566e+00,  5.10620117e+00,  1.02988281e+03,  7.35490723e-01,
#         4.52294922e+02])

# paramValues[13953]=np.array([-2.23205566e+00,  5.10620117e+00,  1.02988281e+03,  7.35490723e-01,
#         4.52294922e+02])
# paramValues[14333]=np.array([-5.11108398e-01,  5.10559082e+01,  2.39628906e+03,  1.94084473e-01,
#         1.05419922e+02])
# paramValues[14334]=np.array([-5.11108398e-01,  1.20653076e+02,  2.04238281e+03,  4.71313477e-02,
#         1.05419922e+02])
# paramValues[15457]=np.array([-6.11938477e-01,  7.73937988e+01,  2.26035156e+03,  5.13127441e-01,
#         2.82470703e+02])

# paramValues[15458]=np.array([-3.32067871e+00,  1.14139404e+02,  2.26035156e+03,  5.13127441e-01,
#         2.82470703e+02])

# paramValues[15459]=np.array([-3.32067871e+00,  7.73937988e+01,  1.97207031e+03,  5.13127441e-01,
#         2.82470703e+02])

# paramValues[15460]=np.array([-3.32067871e+00,  7.73937988e+01,  2.26035156e+03,  3.81643066e-01,
#         2.82470703e+02])
# paramValues[15466] =np.array([-6.11938477e-01,  1.14139404e+02,  1.97207031e+03,  3.81643066e-01,
#         4.82470703e+02])
np.savetxt('Parameters_onset.txt', paramValues, fmt='%f')
print(paramValues.shape)
yArray = np.zeros([paramValues.shape[0],2])
#%%    
for index, param in enumerate(paramValues):
    print(index)
    temp_name='model'+str(index)
    mdir=os.path.join(wdir,temp_name)
    in_dict = {'slope':10**param[0],
               'depth':param[1],
               'width':param[2],
               'aq_thickness':param[4],
               'nx':1000, # use same number of cells per model
               'Lx':5e3,
               'ztype':'spread', # zbot type = ['top','bot','spread']
               'wdir':mdir,
               'rel_level':[param[3]],# not all of them
               'RK_array':[1E-5,1E-4,1E-3,1E-2,1E-1,1][::-1],'max_iter':10
               }
    
    all_RK_list,seep_active_all_RL,data_df = seep_active(**in_dict)
    all_RK_regional_list,seep_active_regonal_all_RL,data_regional_df=seep_active_regional(**in_dict)
    yArray[index]=(seep_active_all_RL[0],seep_active_regonal_all_RL[0])
np.savetxt('outputs_both_onset.txt',yArray,fmt='%f')
#%%
nsamples = np.arange(40, 1001, 20) 
# Arrays to store the index estimates
S1_estimates = np.zeros([problem['num_vars'],len(nsamples)])
ST_estimates = np.zeros([problem['num_vars'],len(nsamples)])
S1r_estimates = np.zeros([problem['num_vars'],len(nsamples)])
STr_estimates = np.zeros([problem['num_vars'],len(nsamples)])
Si = sobol.analyze(problem, yArray[:,1], print_to_console=True)
for i in range(len(nsamples)):
    print('n= '+ str(nsamples[i]))
    # Generate samples
    paramValues=saltelli.sample(problem,nsamples[i])
    # Run model for all samples
    siz=paramValues.shape[0]
    output = yArray[:siz,:]
    # Perform analysis
    results = sobol.analyze(problem, output[:,0],print_to_console=False)
    results_reg =sobol.analyze(problem,output[:,1],print_to_console=False)
    # Store estimates
    ST_estimates[:,i]=results['ST']
    S1_estimates[:,i]=results['S1']
    STr_estimates[:,i]=results_reg['ST']
    S1r_estimates[:,i]=results_reg['S1']
#seep_face_properties(slope=0.1,depth=20,width=1000,rel_lev=rel_level[j],aq_thick=thickness_list[i],r_over_K=0.008,dx=1,nlay=10,wdir=newpath,modelname=temp_name,Lx=8000,rerun=True)
sensArray = np.zeros([8,problem['num_vars']])
#Sobol=Si.to_df
T_Si, first_Si, (idx, second_Si) = sobol.Si_to_pandas_dict(Si)
rk=sobol.analyze(problem, yArray[:,0], print_to_console=True)
rrk = sobol.analyze(problem, yArray[:,1], print_to_console=True)
np.savetxt('Globsens_uzf_RkSt.txt', rk['ST'], fmt='%f')
np.savetxt('Globsens_uzf_RkS1.txt', rk['S1'], fmt='%f')
np.savetxt('Globsens_uzf_RRKSt.txt', rrk['ST'], fmt='%f')
np.savetxt('Globsens_uzf_RRKS1.txt', rrk['S1'], fmt='%f')
#%%
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2)
fig.suptitle("Sobol Indices for Onset of GWD")
fig.set_figheight(26/2.54)
fig.set_figwidth(19/2.54)
ax1.title.set_text("Total Si for (R/K)* ")
ax1.set(xlabel='Number of Samples')
ax1.plot(nsamples,ST_estimates[0,:],label="Slope")
ax1.plot(nsamples,ST_estimates[1,:],label="D/W")
ax1.plot(nsamples,ST_estimates[2,:],label="W/L")
ax1.plot(nsamples,ST_estimates[3,:],label="RL")
ax1.plot(nsamples,ST_estimates[4,:],label="B/L")
ax1.legend()
#%%
ax2.title.set_text("First-order Si for (R/K)* ")
ax2.set(xlabel='Number of Samples')
ax2.plot(nsamples,S1_estimates[0,:],label="Slope")
ax2.plot(nsamples,S1_estimates[1,:],label="D/W")
ax2.plot(nsamples,S1_estimates[2,:],label="W/L")
ax2.plot(nsamples,S1_estimates[3,:],label="RL")
ax2.plot(nsamples,S1_estimates[4,:],label="B/L")
ax2.legend()
#%%
ax3.title.set_text("Total Si for upper (R/K)* ")
ax3.set(xlabel='Number of Samples')
ax3.plot(nsamples,STr_estimates[0,:],label="Slope")
ax3.plot(nsamples,STr_estimates[1,:],label="D/W")
ax3.plot(nsamples,STr_estimates[2,:],label="W/L")
ax3.plot(nsamples,STr_estimates[3,:],label="RL")
ax3.plot(nsamples,STr_estimates[4,:],label="B/L")
ax3.legend()
#%%
ax4.title.set_text("First-order Si for upper (R/K)* ")
ax4.set(xlabel='Number of Samples')
ax4.plot(nsamples,S1r_estimates[0,:],label="Slope")
ax4.plot(nsamples,S1r_estimates[1,:],label="D/W")
ax4.plot(nsamples,S1r_estimates[2,:],label="W/L")
ax4.plot(nsamples,S1r_estimates[3,:],label="RL")
ax4.plot(nsamples,S1r_estimates[4,:],label="B/L")
ax4.legend()