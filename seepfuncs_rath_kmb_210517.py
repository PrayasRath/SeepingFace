

import os
import numpy as np
import pandas as pd
import flopy
from flopy.utils import binaryfile as bf
import flopy.utils.postprocessing as pp
import matplotlib.pyplot as plt
from flopy.plot import PlotCrossSection as pcs

#%%
# need to iterate the model and turn off recharge cells where drns are active and check to see if the next model iteration
# has new drains or not - once the drn cells stabilize, no additional rech to cells need to turn off

def seep_active(slope,depth,width,aq_thickness,nx=None,wdir=None,Lx=None,
                rel_level=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95],
                RK_array=[1E-5,1E-4,1E-3,1E-2,1E-1,1][::-1],
                max_iter=20,ztype='spread'):

 
    seep_active_all_RL=[]
    all_RK_list = []
    all_wt = []
    model_runs=[]
    save_seep_icount = []
    failed_models=[]
    passed_models=[]
    all_fnames = []
    all_data = [ ] # to store as dataframe
    # base dictionary for inputs
    in_dict={'slope':slope,'depth':depth,'width':width,
         'aq_thick':aq_thickness,
         'ztype':ztype,
         'Lx':Lx,'ncol':nx,
         'use_drn':False,'kmult':1}
    icount=-1
    for k in range (len(rel_level)):
        rel_lev_RK_list = []
        all_wt2 = []
        seep_actives=0
        
        # Update inputs
        temp_dict = in_dict.copy()
        temp_dict['rel_lev'] = rel_level[k]
        
        for i in range (len(RK_array)):
            icount += 1
            if (np.log10(RK_array[i])>0):
                temp_name='rel'+str(rel_level[k])+'Rkp'+str(np.log10(RK_array[i]))
            else:
                temp_name ='rel'+str(rel_level[k])+'Rkm'+str(np.abs(np.log10(RK_array[i])))

            wdir2 = os.path.join(wdir,temp_name)

            if not os.path.isdir(wdir2):
                os.makedirs(wdir2)
            
            print("---- Main: {0} ------".format(temp_name))
            
            temp_dict.update({'modelname':temp_name,
                              'wdir':wdir2,
                              'r_over_K':RK_array[i]})
            all_fnames.append(os.path.join(wdir2,temp_name))
            # Run model
            test_result,elev,cell_types = seep_model(**temp_dict)
            # Analyze output
            
            all_out=seep_face_properties(modelname=temp_name,wdir=wdir2,
                                             cell_types=cell_types,
                                             dem=elev,dx=Lx/nx)
            
            dr,dr_seep_total,wt_temp,sff,slf,slf_2d,dr_seep_out_total,ch_ratio,ch_out,q_ratio,sl_regional,sf_regional=all_out
            if (test_result=='failed'):
                failed_models.append(temp_name)
            else:
                passed_models.append(temp_name)
            rel_lev_RK_list.append([icount,rel_level[k],RK_array[i],dr_seep_total])
            all_wt2.append(wt_temp)
            
            # Store more output data
            all_data.append([icount,os.path.join(wdir2,temp_name),slope,depth,width,
                             Lx,nx,aq_thickness,rel_level[k],RK_array[i],ztype,
                             dr_seep_total,sff,slf,dr_seep_out_total,ch_ratio,
                             ch_out,q_ratio,'outer',False]) # assume not onset of seepage for now
            
            if dr_seep_total == 0: # if current model has no seepage from a drain
                seep_actives=RK_array[i-1] # take last model with active seepage
                main_RK = seep_actives
                main_icount=icount
                break
            
            drlast = dr_seep_total
        
        if seep_actives == 0:
            print('No R/K threshold found with values provided.')
            
        
        new_RK_array = np.logspace(np.log10(RK_array[i]),np.log10(RK_array[i-1]),max_iter+2) # includes edges
        new_RK_array = new_RK_array[1:-1]
        new_RK_array = np.round(new_RK_array,decimals=8)
        # dRK = seep_actives/max_iter # can be updated in while loop
        
        for iRK in range(len(new_RK_array)):
            
            icount += 1
            
            seep_actives=new_RK_array[iRK] # needs to come before model run or repeating earlier model
            if (np.log10(seep_actives)>0):
                temp_name='rel'+str(rel_level[k])+'Rkp'+str(np.round(np.log10(seep_actives),decimals=8))
            else:
                temp_name ='rel'+str(rel_level[k])+'Rkm'+str(np.round(np.abs(np.log10(seep_actives)),decimals=8))
            wdir2 = os.path.join(wdir,temp_name)
            if not os.path.isdir(wdir2):
                os.makedirs(wdir2)
            print("---- Inner: {0} ------".format(temp_name))
            all_fnames.append(os.path.join(wdir2,temp_name))
            temp_dict.update({'modelname':temp_name,
                              'wdir':wdir2,
                              'r_over_K':seep_actives})
            
            # Run model
            test_result,elev,cell_types = seep_model(**temp_dict)
            # Analyze output
            
            all_out=seep_face_properties(modelname=temp_name,wdir=wdir2,
                                             cell_types=cell_types,
                                             dem=elev,dx=Lx/nx)
            
            dr,dr_seep_total,wt_temp,sff,slf,slf_2d,dr_seep_out_total,ch_ratio,ch_out,q_ratio,sl_regional,sf_regional=all_out
            
            
            rel_lev_RK_list.append([icount,rel_level[k],new_RK_array[iRK],dr_seep_total])
            all_wt2.append(wt_temp)
            if test_result=='failed':
                    failed_models.append(temp_name)
            else:
                passed_models.append(temp_name)
            # print(rel_level[k],seep_actives)
            
            
            
            if dr_seep_total > 0.:
                save_RK = seep_actives
                save_seep_icount.append(icount)
                all_data.append([icount,os.path.join(wdir2,temp_name),slope,depth,width,
                             Lx,nx,aq_thickness,rel_level[k],new_RK_array[iRK],ztype,
                             dr_seep_total,sff,slf,dr_seep_out_total,ch_ratio,
                             ch_out,q_ratio,'inner',True])
                break # no need to keep running models for this rel_level scenario
            elif iRK == len(new_RK_array)-1:
                save_RK = main_RK # could be in between last value in new_RK_array and RK_array[i-1]
                save_seep_icount.append(main_icount)
                all_data.append([icount,os.path.join(wdir2,temp_name),slope,depth,width,
                             Lx,nx,aq_thickness,rel_level[k],new_RK_array[iRK],ztype,
                             dr_seep_total,sff,slf,dr_seep_out_total,ch_ratio,
                             ch_out,q_ratio,'inner',False])
            else:
                all_data.append([icount,os.path.join(wdir2,temp_name),slope,depth,width,
                             Lx,nx,aq_thickness,rel_level[k],new_RK_array[iRK],ztype,
                             dr_seep_total,sff,slf,dr_seep_out_total,ch_ratio,
                             ch_out,q_ratio,'inner',False])
            
            
            
            # drlast = dr_seep_total
            
            # if icount > max_iter:
            #     print("----- Max iterations exceeded ------")
            
            
        print('Active drains at R/K={0:3.2e}'.format(save_RK))
        all_RK_list.append(np.array(rel_lev_RK_list))
        seep_active_all_RL.append(save_RK) # activation of seepage at this R/K value
        all_wt.append(np.array(all_wt2))
        model_runs.append(icount) # isn't this just a list of higher and higher numbers in order? why not just return icount?
    
    cols = ['model_num','fpath','regional_slope','valley_depth','valley_width',
            'domain_L','nx','aq_thick_b','rel_level','roverK','zbot_type',
            'valley_seep_Q','SFF','SLF','all_seep_Q','ch_ratio','ch_Q',
            'Q_ratio','model_loop_type','onset_flag']
    all_data_df = pd.DataFrame(all_data,columns=cols)
    
    # Update loops that were the onset of seepage
    all_data_df.loc[all_data_df['model_num'].isin(save_seep_icount),'onset_flag']=True
    
    # need to return these to make them available for later
    return all_RK_list, seep_active_all_RL,all_data_df
#%%
def seep_active_regional(slope,depth,width,aq_thickness,nx=None,wdir=None,Lx=None,
                rel_level=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95],
                RK_array=[1E-5,1E-4,1E-3,1E-2,1E-1,1][::-1],
                max_iter=20,ztype='spread'):

 
    seep_active_all_RL=[]
    all_RK_list = []
    all_wt = []
    model_runs=[]
    save_seep_icount = []
    failed_models=[]
    passed_models=[]
    all_fnames = []
    all_data = [ ] # to store as dataframe
    # base dictionary for inputs
    in_dict={'slope':slope,'depth':depth,'width':width,
         'aq_thick':aq_thickness,
         'ztype':ztype,
         'Lx':Lx,'ncol':nx,
         'use_drn':False,'kmult':1}
    icount=-1
    for k in range (len(rel_level)):
        rel_lev_RK_list = []
        all_wt2 = []
        seep_actives=0
        
        # Update inputs
        temp_dict = in_dict.copy()
        temp_dict['rel_lev'] = rel_level[k]
        
        for i in range (len(RK_array)):
            icount += 1
            if (np.log10(RK_array[i])>0):
                temp_name='rel'+str(rel_level[k])+'Rkp'+str(np.log10(RK_array[i]))
            else:
                temp_name ='rel'+str(rel_level[k])+'Rkm'+str(np.abs(np.log10(RK_array[i])))

            wdir2 = os.path.join(wdir,temp_name)

            if not os.path.isdir(wdir2):
                os.makedirs(wdir2)
            
            print("---- Main: {0} ------".format(temp_name))
            
            temp_dict.update({'modelname':temp_name,
                              'wdir':wdir2,
                              'r_over_K':RK_array[i]})
            all_fnames.append(os.path.join(wdir2,temp_name))
            # Run model
            test_result,elev,cell_types = seep_model(**temp_dict)
            # Analyze output
            
            all_out=seep_face_properties(modelname=temp_name,wdir=wdir2,
                                             cell_types=cell_types,
                                             dem=elev,dx=Lx/nx)
            
            dr,dr_seep_total,wt_temp,sff,slf,slf_2d,dr_seep_out_total,ch_ratio,ch_out,q_ratio,sl_regional,sf_regional=all_out
            if (test_result=='failed'):
                failed_models.append(temp_name)
            else:
                passed_models.append(temp_name)
            rel_lev_RK_list.append([icount,rel_level[k],RK_array[i],dr_seep_out_total])
            all_wt2.append(wt_temp)
            
            # Store more output data
            all_data.append([icount,os.path.join(wdir2,temp_name),slope,depth,width,
                             Lx,nx,aq_thickness,rel_level[k],RK_array[i],ztype,
                             dr_seep_total,sff,slf,dr_seep_out_total,ch_ratio,
                             ch_out,q_ratio,'outer',False]) # assume not onset of seepage for now
            
            if dr_seep_out_total == 0: # if current model has no seepage from a drain
                seep_actives=RK_array[i-1] # take last model with active seepage
                main_RK = seep_actives
                main_icount=icount
                break
            
            drlast = dr_seep_out_total
        
        if seep_actives == 0:
            print('No R/K threshold found with values provided.')
            
        
        new_RK_array = np.logspace(np.log10(RK_array[i]),np.log10(RK_array[i-1]),max_iter+2) # includes edges
        new_RK_array = new_RK_array[1:-1]
        new_RK_array = np.round(new_RK_array,decimals=8)
        # dRK = seep_actives/max_iter # can be updated in while loop
        
        for iRK in range(len(new_RK_array)):
            
            icount += 1
            
            seep_actives=new_RK_array[iRK] # needs to come before model run or repeating earlier model
            if (np.log10(seep_actives)>0):
                temp_name='rel'+str(rel_level[k])+'Rkp'+str(np.round(np.log10(seep_actives),decimals=8))
            else:
                temp_name ='rel'+str(rel_level[k])+'Rkm'+str(np.round(np.abs(np.log10(seep_actives)),decimals=8))
            wdir2 = os.path.join(wdir,temp_name)
            if not os.path.isdir(wdir2):
                os.makedirs(wdir2)
            print("---- Inner: {0} ------".format(temp_name))
            all_fnames.append(os.path.join(wdir2,temp_name))
            temp_dict.update({'modelname':temp_name,
                              'wdir':wdir2,
                              'r_over_K':seep_actives})
            
            # Run model
            test_result,elev,cell_types = seep_model(**temp_dict)
            # Analyze output
            
            all_out=seep_face_properties(modelname=temp_name,wdir=wdir2,
                                             cell_types=cell_types,
                                             dem=elev,dx=Lx/nx)
            
            dr,dr_seep_total,wt_temp,sff,slf,slf_2d,dr_seep_out_total,ch_ratio,ch_out,q_ratio,sl_regional,sf_regional=all_out
            
            
            rel_lev_RK_list.append([icount,rel_level[k],new_RK_array[iRK],dr_seep_total])
            all_wt2.append(wt_temp)
            if test_result=='failed':
                    failed_models.append(temp_name)
            else:
                passed_models.append(temp_name)
            # print(rel_level[k],seep_actives)
            
            
            
            if dr_seep_out_total > 0.:
                save_RK = seep_actives
                save_seep_icount.append(icount)
                all_data.append([icount,os.path.join(wdir2,temp_name),slope,depth,width,
                             Lx,nx,aq_thickness,rel_level[k],new_RK_array[iRK],ztype,
                             dr_seep_total,sff,slf,dr_seep_out_total,ch_ratio,
                             ch_out,q_ratio,'inner',True])
                break # no need to keep running models for this rel_level scenario
            elif iRK == len(new_RK_array)-1:
                save_RK = main_RK # could be in between last value in new_RK_array and RK_array[i-1]
                save_seep_icount.append(main_icount)
                all_data.append([icount,os.path.join(wdir2,temp_name),slope,depth,width,
                             Lx,nx,aq_thickness,rel_level[k],new_RK_array[iRK],ztype,
                             dr_seep_total,sff,slf,dr_seep_out_total,ch_ratio,
                             ch_out,q_ratio,'inner',False])
            else:
                all_data.append([icount,os.path.join(wdir2,temp_name),slope,depth,width,
                             Lx,nx,aq_thickness,rel_level[k],new_RK_array[iRK],ztype,
                             dr_seep_total,sff,slf,dr_seep_out_total,ch_ratio,
                             ch_out,q_ratio,'inner',False])
            
            
            
            # drlast = dr_seep_total
            
            # if icount > max_iter:
            #     print("----- Max iterations exceeded ------")
            
            
        print('Active drains at R/K={0:3.2e}'.format(save_RK))
        all_RK_list.append(np.array(rel_lev_RK_list))
        seep_active_all_RL.append(save_RK) # activation of seepage at this R/K value
        all_wt.append(np.array(all_wt2))
        model_runs.append(icount) # isn't this just a list of higher and higher numbers in order? why not just return icount?
    
    cols = ['model_num','fpath','regional_slope','valley_depth','valley_width',
            'domain_L','nx','aq_thick_b','rel_level','roverK','zbot_type',
            'valley_seep_Q','SFF','SLF','all_seep_Q','ch_ratio','ch_Q',
            'Q_ratio','model_loop_type','onset_flag']
    all_data_df = pd.DataFrame(all_data,columns=cols)
    
    # Update loops that were the onset of seepage
    all_data_df.loc[all_data_df['model_num'].isin(save_seep_icount),'onset_flag']=True
    
    # need to return these to make them available for later
    return all_RK_list, seep_active_all_RL,all_data_df
def make_zbot(dem,aq_thick,ztype='top_thick',nlay=10,min_thick=10):
    
    nrow,ncol = dem[None,:].shape
    zbot=np.zeros((nlay,nrow,ncol)) 
    if ztype in ['top','top_thick']:
        # Original zbot approach, layers below valley bottom
        thickness=(aq_thick-min_thick)/nlay
        zbot[0]=min(dem)-min_thick
        for i in range (1,nlay):
            zbot[i,:,:]=zbot[0]-thickness*i
            
    elif ztype in ['bot','bot_thick']:
        # equal layer thickness from land surface, deepest layer thick
        
        # multiple layers above b
        lay_thick = aq_thick/nlay
         
        for i in range(nlay):
            zbot[i,:,:]=dem-(lay_thick*(i+1))
        
        zbot[-1] = zbot[-1,0,0] # flat bottom set by valley elevation and aq_thick
    elif ztype in ['spread']:
        # each column has layer thickness set by topography
        min_z = dem[0]-aq_thick
        aq_thick_array = (dem-min_z)/nlay # thickness per column
        zbot_depth = np.arange(1,nlay+1)[:,None] * aq_thick_array[None,:]    
        zbot = dem-zbot_depth
        zbot = zbot[:,None,:]
    return zbot

def seep_model(slope=None, depth=None, width=None,rel_lev=None,
                aq_thick=None,r_over_K=None,wdir=None,
                K=1.0, ss=1E-5,sy=0.15,ncol=1000, Lx=1e4, Ly=1e0,
                 nlay=10,nrow=1,z0=0,modelname='2d_res', rerun=False,
                 drn_bool=None,min_thick=1,ztype='top',use_drn=False,kmult=1):
    
    dx = Lx/ncol
    x=np.arange(0,Lx,dx)
    width_act=width/2.
    width_h=int(width_act/(dx))

    res_loc = 0 #int(Lx/dx/2.)-1
    # lake_loc=int(res_loc/dx)-1 #domain world
    dem1=np.zeros(x.shape[0],dtype=np.float64)
    
    # Make reservoir depression
    dem1[res_loc:res_loc+width_h+1] = np.linspace(z0-depth,z0,width_h+1)
    
    # Make upslope topography
    dem1[res_loc+width_h+1:] = (x[res_loc+width_h+1:]-x[res_loc+width_h])*slope + z0
    
    
    # dem=dem1.reshape(1,dem1.shape[0])
    
    lake_head = rel_lev*depth + dem1[res_loc]

    # 

    delc=Ly/nrow
    delr=Lx/ncol
    
    reservoir_bound=np.zeros_like(x).astype(bool)
    # i=0
    j=0
    
    for j in range(width_h +2):
        if (dem1[res_loc+j]>lake_head):
        
           break
        right_bound=res_loc+j

 

    reservoir_bound[:right_bound+1]=True
    
    cell_types = np.ones_like(x) # 1 = normal active cell
    cell_types[res_loc:res_loc+width_h+1] = 2 # 2 = waterbody valley
    cell_types[reservoir_bound] = -2 # -2 = surface waterbody
    
    cbc_fname = os.path.join(wdir,modelname+'.cbc')        
    if not os.path.isfile(cbc_fname) or rerun: # only run model if it doesn't already have outputs
        mf = flopy.modflow.Modflow(modelname, exe_name='MODFLOW-NWT_64',version='mfnwt',model_ws=wdir)
    
        
        zbot = make_zbot(dem1,aq_thick,ztype=ztype,nlay=nlay,min_thick=min_thick)
  
        nper=1
        perlen=1
        steady=[True]
        nstp=1
        # Use default for units: meters and days
        dis=flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, 
                                     delc=delc,top=dem1[None,:],nper=nper,perlen=perlen,
                                     nstp=nstp, botm=zbot, steady=steady)

        # BAS package        
        ibound=np.ones((nlay,nrow, ncol))
        ibound[0,:,cell_types==-2]=-1 # waterbody as constant head

        ihead = lake_head * np.ones((nlay, nrow, ncol), np.float)
        ihead[0,:,cell_types==-2]=lake_head # waterbody constant head value

        bas= flopy.modflow.ModflowBas(mf, ibound=ibound, strt=ihead)
    
        upw = flopy.modflow.ModflowUpw(mf,laytyp=1,hk=K,ss=ss,sy=sy,vka=K,ipakcb=53)
        
        drain=np.where(cell_types!=-2)[0]  # everywhere but the constant heads
    
        dr_col=drain.astype(int)
        dr_row=[0]*len(dr_col)
        dr_layer=[0]*len(dr_col)
        # z_down=zbot[0][0][0]
        # dz=dem1[reservoir_bound==False]-z_down
#        cond3=[(1e3*K*delr*delc)/(thickness/2)]*len(drx)
        top_lay_thick = dem1[dr_col]-zbot[0,0,dr_col]
        cond3=(K*delr*delc)/(top_lay_thick/2) # vertical conductance from cell center to drain
        elev_vals=dem1[drain].copy()
        idrn=np.column_stack([dr_layer,dr_row,dr_col,elev_vals,cond3])
        
        #recharge per year
        rech=K*r_over_K # defined by ratio
        rech_array = rech*np.ones_like(cell_types,dtype=float)
        
        if drn_bool is not None:
            rech_array[drn_bool | (cell_types==-2)] = 0 # turn off recharge for cells with active drainage or constant heads 
        
        if use_drn:
            
            drn=flopy.modflow.ModflowDrn(mf,stress_period_data=idrn,ipakcb=53,options=['NOPRINT']) 
            rch = flopy.modflow.ModflowRch(mf,nrchop=3,rech=rech_array[None,:],ipakcb=53)
        else:
            iuzfbnd = np.ones((nrow,ncol),dtype=int)
            iuzfbnd[0,cell_types==-2] = 0
            K_array = K*np.ones((nrow,ncol),dtype=float)
            K_array[:,cell_types==1] = K*kmult # outside of valley
            uzf = flopy.modflow.ModflowUzf1(mf,nuztop=1,surfdep=0.01*dx,iuzfopt=1,
                                            vks=K_array,finf=rech_array[None,:],iuzfbnd=iuzfbnd,ipakcb=53,
                                            )
    
        spd = {(0, 0): [ 'save head', 'save budget']}
        oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)
        # Add nwt package to the MODFLOW model
        # pcg= flopy.modflow.ModflowPcg(mf)
        nwt = flopy.modflow.ModflowNwt(mf,headtol=0.0001,fluxtol=500,options='COMPLEX',maxiterout=2000,iprnwt=1)
    
        # Add Modpath to the model
        
        # Write the MODFLOW model input files
        mf.write_input()
    
        # Run the MODFLOW model
        success, buff = mf.run_model(silent=True,report=True)
        if buff[-1]=='  Normal termination of simulation':
                test_result='pass'
        else:
            test_result='failed'
        
        return test_result,dem1,cell_types
    else:
        return 'exists',dem1,cell_types

def plot_model(modelname=None,wdir=None,ax=None):
    
    mf = flopy.modflow.Modflow.load('{}.nam'.format(modelname),model_ws=wdir,exe_name='mfnwt')
    xsec = pcs(mf,line={'row':0},ax=ax)
    xsec.plot_grid(colors='k')
    # xsec.plot_bc('DRN',color='r')
    
    wt,h=load_wt(modelname,wdir)
    csa = xsec.plot_array(h,head=h,masked_values=[-1e+30])
    xsec.plot_surface(wt[None,:],color='b')

    xsec.plot_ibound()
    cb = plt.colorbar(csa, shrink=0.75)

def seep_test(modelname=None,wdir=None,cell_types=None):
    
    cbc_fname = os.path.join(wdir,'{}.cbc'.format(modelname))
    cbb = bf.CellBudgetFile(cbc_fname)
    dr=cbb.get_data(text='DRAINS',totim=1)
    cbb.close()
    cbb = None
    
    active_drains = dr[0]['node'][dr[0]['q']<0] - 1 # column index
    
    drn_bool = np.zeros_like(cell_types,dtype=bool)
    drn_bool[active_drains] = True
    
    return drn_bool

def load_wt(modelname=None,wdir=None,totim=1):
    hds = bf.HeadFile(os.path.join(wdir,'{}.hds'.format(modelname)))
    head = hds.get_data(totim=totim)
    water_table=pp.get_water_table(head,-1e+30)
    hds.close()
    hds = None
    return water_table,head

def seep_face_properties(modelname=None,wdir=None,
                         cell_types=None,use_drn=False,totim=None,
                         dem=None,dx=None,plot_bool=False,
                         seep_depth_threshold=-1E-2):
    
    cbc_fname = os.path.join(wdir,'{}.cbc'.format(modelname))
    cbb = bf.CellBudgetFile(cbc_fname)
    ch=cbb.get_data(text='CONSTANT HEAD',totim=totim)
    ch_inds = ch[0]['node']-1 # back to 0 base
    
    res_water = ch_inds.copy()
    res_valley = (np.abs(cell_types)==2).nonzero()[0]
    max_seep = np.sort(np.asarray(list(set(res_valley)-set(res_water)))) # indexes for maximum seepage extent - unflooded valley inds
    
    
    
    if use_drn:
        dr=cbb.get_data(text='DRAINS',totim=totim)
        drain_inds = dr[0]['node']-1 # back to 0 base
        active_drains = dr[0]['node'][dr[0]['q']<0] - 1 # column index
        
        dr_seep = np.in1d(dr[0]['node'],max_seep+1) # locate drain assignments in potential seepage zone
        dr_seep_total = sum(-dr[0]['q'][dr_seep]) # total groundwater discharge from valley seeps
        
        # Regional drain discharge
        dr_seep_outside = ~dr_seep # locate drains outside of waterbody valley
        dr_seep_out_total = sum(-dr[0]['q'][dr_seep_outside])
        
        rch_in = cbb.get_data(text='RECHARGE',totim=totim)
        r_inds = (rch_in[0][1].squeeze()>0).nonzero()[0]
    else:
        # Load from uzf
        
        rch_in = cbb.get_data(text='UZF RECHARGE')[0].squeeze().sum(axis=0) # add rech to other layers
        dr = cbb.get_data(text='SURFACE LEAKAGE')[0].squeeze()
        dr_seep = (dr<0) & (cell_types==2) # where seepage occurs in valley, boolean
        dr_seep_total = np.sum(-dr[dr_seep])
        
        # Regional drain discharge
        dr_seep_outside = (dr<0) & (cell_types==1) # locate drains outside of waterbody valley
        dr_seep_out_total = np.sum(-dr[dr_seep_outside])
    
    cbb.close()
    cbb = None
    
    water_table, head=load_wt(modelname,wdir)
    
    
    
    # Constant head groundwater flow
    ch_seep_total=sum(-ch[0]['q'])
    chs=-ch[0]['q'] # constant head flow positive into surface water, negative for flow into aquifer
    ch_in=sum(chs[chs>0]) # total baseflow to surface water
    ch_out=sum(chs[chs<0]) # total recharge from surface water
    
    # Waterbody valley drain discharge
    dr_seep_len = np.sum(dr_seep)*dx # active seepage face length 
    
    
    
    
    # dr_ch_ratio=dr_seep_total/ch_seep_total
    
    valley_bound=np.zeros_like(cell_types).astype(bool)
    valley_bound[res_valley]=True
    
    seep_lngth=np.squeeze(water_table-dem)>=seep_depth_threshold  #less than head tol
    seep_length_val_bool=seep_lngth[valley_bound]
    seep_length_val=(len(seep_length_val_bool.nonzero()[0])-len(res_water))*dx
    seep_length_all=(len(seep_lngth.nonzero()[0])-len(res_water))*dx
    slf=seep_length_val/((len(res_valley)-len(res_water))*dx)
    slf_regional=(seep_length_all-seep_length_val)/(dx*(len(dem)-len(res_valley)))
#    max_seep=np.asarray(list(set(res_valley)-set(res_water)))
    seep_lngth_bool=seep_lngth[list(set(res_valley)-set(res_water))]
    actual_seep=max_seep[seep_lngth_bool]
    total_seep_len=0
    lower_bank_seep=0
    upper_bank_seep=0
    if len(actual_seep)!=0:
         dif=np.diff(actual_seep)
         if len(set(dif))>1:
            diff=np.where(dif!=1)[0][0]
            lower_bank_seep=0
            upper_bank_seep=0
            sq_dist=0
            for i in range(actual_seep[0],actual_seep[diff]+1):
                sq_dist=sq_dist+((((dem[i+1]-dem[i])**2+(dx)**2))**0.5)
            lower_bank_seep=sq_dist
            sq_dist=0
            for i in range(actual_seep[diff+1],actual_seep[-1]):
                sq_dist=sq_dist+((((dem[i+1]-dem[i])**2+(dx)**2))**0.5)
            upper_bank_seep=sq_dist
            total_seep_len=upper_bank_seep+lower_bank_seep
         else:
             total_seep_len=0 
             sq_dist=0
             for i in range(actual_seep[0],actual_seep[-1]+1):
                 sq_dist=sq_dist+((((dem[i+1]-dem[i])**2+(dx)**2))**0.5)
             total_seep_len=sq_dist 
    max_seep_len_l=0
    max_seep_len_r=0
    #    dif=np.diff(max_seep)
    #    diff=np.where(dif!=1)[0][0]
    sq_dist=0
    for i in range (max_seep[0],max_seep[-1]):
                sq_dist=sq_dist+((((dem[i+1]-dem[i])**2+(dx)**2))**0.5)
    max_seep_len=sq_dist
    #    sq_dist=0
    #    for i in range (max_seep[diff+1],max_seep[-1]):           
    #                sq_dist=sq_dist+((((dem1[i+1]-dem1[i])**2+(x[i+1]-x[i])**2))**0.5)
            #    plt.step(x,water_table, where='post')
    #    plt.step(x,dem1, where='post')
    if plot_bool:
        x = np.arange(len(cell_types))*dx
        plt.plot(x,water_table,label=modelname)
        plt.plot(x,dem)
        plt.legend()
        plt.show()
    #    max_seep_len_r=sq_dist
    #    max_seep_len=max_seep_len_r+max_seep_len_l
    #plt.plot(x,dem1-water_table)
    # print(total_seep_len,dr_ch_ratio,ch_seep_total,dr_seep_total)
    # return slf,seep_length_val,total_seep_len,dr_ch_ratio,ch_seep_total,dr_seep_total, dr_seep_out_total
    sff=dr_seep_total/(dr_seep_total+dr_seep_out_total+ch_in)
    sff_regional=dr_seep_out_total/((dr_seep_total+dr_seep_out_total+ch_in))
    if max_seep_len > 0:
        slf_2d=total_seep_len/max_seep_len
    else:
        slf_2d = 0
        
    if ch_in==0:
        ch_ratio=-8888
        q_ratio=-8888
    else:    
        ch_ratio=ch_out/ch_in
    q_ratio=ch_out/(dr_seep_total+ch_in)
    #print(dr_seep_total,ch_seep_total,sff,slf,slf_2d,dr_seep_out_total,total_seep_len,max_seep_len)
    #print(sff,slf,slf_2d,total_seep_len,max_seep_len,dr_seep,dr_seep_len,len(res_valley),len(res_water),actual_seep,max_seep)
    return -dr[0],dr_seep_total,water_table,sff,slf,slf_2d,dr_seep_out_total,ch_ratio,ch_out,q_ratio,sff_regional,slf_regional
    
    
