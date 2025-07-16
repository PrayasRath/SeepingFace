# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 07:36:11 2021

@author: prath1
"""
import os
import numpy as np
# from skimage import io
import flopy
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
# import flopy.utils.geometry as geo
import flopy.utils.postprocessing as pp
import math
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#
#os.chdir("E:\\August_20\\Seepage_onset")
wdir = os.path.join(r'E:\seepage_face_UZF2','SLF_SFF')
os.chdir(wdir)
#%%
#df5=pd.read_csv("drainage_basedseepage_avtivation_10layer_diff_gwd2.csv",header=None)
df1=pd.read_csv("uzf_diff_RK_SLF_low_s.csv",header=None)
df2=pd.read_csv("uzf_diff_RK_SLF_regional_low_s.csv",header=None)
#%%
df3=pd.read_csv("uzf_diff_slope_SLF_mid_r.csv",header=None)
df4=pd.read_csv("uzf_diff_slope_SLF_regional_mid_r.csv",header=None)
#%%
df5=pd.read_csv("uzf_diff_vslope_SLF_low_s.csv",header=None)
df6=pd.read_csv("uzf_diff_vslope_SLF_regional_low_s.csv",header=None)
#%%
df7=pd.read_csv("uzf_diff_width_SLF_low_s.csv",header=None)
df8=pd.read_csv("uzf_diff_width_SLF_regional_low_s.csv",header=None)
#%%
df9=pd.read_csv("uzf_diff_thick_SLF_low_s.csv",header=None)
df10=pd.read_csv("uzf_diff_thick_SLF_regional_low_s.csv",header=None)
#%%
from scipy.interpolate import griddata
#fig, ((ax1, ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(nrows=2,ncols=3)
fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8),(ax9,ax10)) = plt.subplots(nrows=5,ncols=2)
fig.suptitle("Seepage Flow Fraction For S=0.001")
fig.set_figheight(12)
fig.set_figwidth(8)
#%%
X=[10,20,30,40,50,60,70,80,90]
Y=np.log10(np.array([1e-4,10**(-3.5),1e-3,10**(-2.5),1e-2,10**-(1.5),1e-1,10**(-0.5),10**(0)]))

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()

#ZZ=griddata(points,np.log10(df2.values[:,0:11]),(XX,YY),method='linear')
im=ax1.contourf(X,Y,(df1.values),levels=100,vmin=0,vmax=1,cmap='viridis_r')
ax1.title.set_text('SLFL')
scat=ax1.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax1.contour(XX,YY,(df1.values),np.arange(0.1,1,0.1),colors='black')
ax1.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax1.set(ylabel="log10(R/K)) [-]")

#ax1.set(xlabel="Relative levels [%]")
#%%
X=[10,20,30,40,50,60,70,80,90]
Y=np.log10(np.array([1e-4,10**(-3.5),1e-3,10**(-2.5),1e-2,10**-(1.5),1e-1,10**(-0.5),10**(0)]))

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()

#ZZ=griddata(points,np.log10(df2.values[:,0:11]),(XX,YY),method='linear')
im=ax2.contourf(X,Y,(df2.values),levels=100,vmin=0,vmax=1,cmap='viridis_r')
ax2.title.set_text('SLFU')
scat=ax2.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax2.contour(XX,YY,(df2.values),np.arange(0.1,1,0.1),colors='black')
ax2.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax2.set(ylabel="log10(R/K) [-]")
#ax2.set_aspect('equal')
#ax2.set(xlabel="Relative levels [%]")
#%%

X=[10,20,30,40,50,60,70,80,90]
Y=np.log10(np.array([10**(-5),10**(-4.5),10**(-4),10**(-3.5),10**(-3),10**(-2.5),10**(-2),10**(-1.5),10**(-1),10**(-0.5)]))

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()

#ZZ=griddata(points,np.log10(df2.values[:,0:11]),(XX,YY),method='linear')
im=ax3.contourf(X,Y,(df3.values),levels=100,vmin=0,vmax=1,cmap='viridis_r')
#ax3.title.set_text('Effect of upper slope')
scat=ax3.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax3.contour(XX,YY,(df3.values),np.arange(0.1,1,0.1),colors='black')
ax3.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax3.set(ylabel="log10(S) [-]")
#ax3.set_aspect('equal')
#ax3.set(xlabel="Relative levels [%]")
#%%
X=[10,20,30,40,50,60,70,80,90]
Y=np.log10(np.array([10**(-5),10**(-4.5),10**(-4),10**(-3.5),10**(-3),10**(-2.5),10**(-2),10**(-1.5),10**(-1),10**(-0.5)]))

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()

#ZZ=griddata(points,np.log10(df2.values[:,0:11]),(XX,YY),method='linear')
im=ax4.contourf(X,Y,df4,levels=100,vmin=0,vmax=1,cmap='viridis_r')
#ax4.title.set_text('Diffrence in order of magnitude for various upper slopes')
scat=ax4.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax4.contour(XX,YY,df4,np.arange(0.1,1,0.1),colors='black')
ax4.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax4.set(ylabel="log10(S) [-]")
#ax4.set_aspect('equal')
#ax4.set(xlabel="Relative levels [%]")
#%%
X=[10,20,30,40,50,60,70,80,90]
Y=(np.array([0.01,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]))

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()
#Z=df9.values.flatten()
#ZZ=griddata(points,np.log10(df2.values[:,0:11]),(XX,YY),method='linear')
im=ax5.contourf(X,Y,(df5.values),levels=100,vmin=0,vmax=1,cmap='viridis_r')
#ax5.title.set_text(' Effect of lower slope')
scat=ax5.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax5.contour(XX,YY,df5.values,np.arange(0.1,1,0.1),colors='black')
ax5.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax5.set(ylabel="d/w [-]")
#ax5.set_aspect('equal')
#ax5.set(xlabel="Relative levels [%]")
#%% width_low_slope
X=[10,20,30,40,50,60,70,80,90]
Y=np.array([0.01,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3])

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()
#Z=df9.values.flatten()
#ZZ=griddata(points,np.log10(df2.values[:,0:11]),(XX,YY),method='linear')
im=ax6.contourf(X,Y,(df6.values),levels=100,vmin=0,vmax=1,cmap='viridis_r')
#ax6.title.set_text('Diffrence in order of magnitude for various lower slopes')
scat=ax6.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax6.contour(XX,YY,df6.values,np.arange(0.1,1,0.1),colors='black')
ax6.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax6.set(ylabel="d/w [-]")
#ax6.set_aspect('equal')
#%% width_low_slope
X=[10,20,30,40,50,60,70,80,90]
Y=np.array([0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5])

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()
#Z=df5.values.flatten()
#ZZ=griddata(points,np.log10(Z),(XX,YY),method='linear')
im=ax7.contourf(X,Y,df7.values,levels=100,vmin=0,vmax=1,cmap='viridis_r')
#ax7.title.set_text('Effect of width')
scat=ax7.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax7.contour(XX,YY,df7.values,np.arange(0.1,1,0.1),colors='black')
ax7.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax7.set(ylabel="w/L [-]")
#ax7.set_aspect('equal')
#%% width_low_slope
X=[10,20,30,40,50,60,70,80,90]
Y=np.array([0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5])

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()
#Z=df9.values.flatten()
#ZZ=griddata(points,np.log10(df2.values[:,0:11]),(XX,YY),method='linear')
im=ax8.contourf(X,Y,df8.values,levels=100,vmin=0,vmax=1,cmap='viridis_r')
#ax8.title.set_text('Log(Regional Onset)-Log(Valley Onset)')
scat=ax8.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax8.contour(XX,YY,df8.values,np.arange(0.1,1,0.1),colors='black')
ax8.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax8.set(ylabel="w/L [-]")
#ax8.set_aspect('equal')
#%% width_low_slope
X=[10,20,30,40,50,60,70,80,90]
Y=[0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15]

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()
#Z=df9.values.flatten()
#ZZ=griddata(points,np.log10(df2.values[:,0:11]),(XX,YY),method='linear')
im=ax9.contourf(X,Y,df9.values,levels=100,vmin=0,vmax=1,cmap='viridis_r')
#ax9.title.set_text('Effect of aquifer thickness')
scat=ax9.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax9.contour(XX,YY,df9.values,np.arange(0.1,1,0.1),colors='black')
ax9.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax9.set(ylabel="B/L [-]",xlabel='Relative levels[%]')
#ax9.set_aspect('equal')
#%% width_low_slope
X=[10,20,30,40,50,60,70,80,90]
Y=[0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15]

XX,YY=np.meshgrid(X,Y)
length_values=len(X)*len(Y)
points = np.empty((length_values, 2))
points[:, 0] = XX.flatten()
points[:, 1] = YY.flatten()
#Z=df9.values.flatten()
#ZZ=griddata(points,np.log10(df2.values[:,0:11]),(XX,YY),method='linear')
im=ax10.contourf(X,Y,df10.values,levels=100,vmin=0,vmax=1,cmap='viridis_r')
#ax10.title.set_text('Effect of aquifer thickness')
scat=ax10.scatter(XX,YY,facecolors='none', edgecolors='black',s=10)
cl=ax10.contour(XX,YY,df10.values,np.arange(0.1,1,0.1),colors='black')
ax10.clabel(cl,cl.levels,fontsize=10,inline=False,fmt='%1.2f',colors='w')
ax10.set(ylabel="B/L [-]",xlabel='Relative levels[%]')
#ax10.set_aspect('equal')
#%%
fig.subplots_adjust(right=0.9)
#cbar=fig.colorbar(im,ax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9],orientation='horizontal')
import matplotlib as mpl
cmap=mpl.cm.viridis_r
norm= mpl.colors.Normalize(vmin=0, vmax=1)
cbar_ax = fig.add_axes([0.92, 0.15, 0.025, 0.7])
cbar=mpl.colorbar.ColorbarBase(cbar_ax,cmap=cmap,norm=norm)#orientation='horizontal')
#cbar=fig.colorbar(im,orientation='horizontal')
#cbar.set_ticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
cbar.set_label("SLF[-]")
plt.show()