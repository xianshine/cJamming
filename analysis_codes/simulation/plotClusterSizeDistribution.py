#!/usr/bin/python

import re
import sys
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
#mpl.use('GTK')
import matplotlib.pyplot as plt
import math as m

from scipy.interpolate import spline
from scipy.interpolate import BSpline
import scipy.interpolate as si

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar


#plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['svg.fonttype'] = 'svgfont'

# set tick width
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 1

mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 1

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'

mpl.rcParams['xtick.major.pad']='6'
mpl.rcParams['ytick.major.pad']='6'
mpl.rcParams['axes.labelpad']='2'
#mpl.rcParams['axes.titlepad']='2'

border=np.fromfile('border2.dat',sep=" ").reshape(35,2)
fiCin=border[:,0]
nCCin=border[:,1]
print border

datRoot='../data/HS2/'
datType=np.dtype([('n',int)])
xn=np.arange(25)-0.5
xp=np.arange(24)

lSs=['b', 'g', 'r', 'c', 'm', 'y','C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']

#nCs = np.arange(20,40,5)
nCs = np.array([20,35])
nCs = np.append(nCs,np.arange(40,120,10))

#rhos=[0.57] ; nCs = np.arange(35,45)
#rhos=[0.63] ; nCs = np.arange(35,36,2) ; nCs=np.append(nCs,np.array([50,60,70,80,90,100]))
#rhos=[0.45] ; nCs = np.arange(20,31)
#rhos=[0.50] ; nCs = np.arange(25,36)
#rhos=[0.55] ; nCs = np.arange(30,41)

rhos=np.arange(0.35,0.691,0.01)
rhos=np.arange(0.45,0.641,0.01)
rhos=np.array([0.45,0.50,0.55,0.57,0.60,0.64])
#rhos=np.array([0.60])

#rhos=np.arange(0.55,0.601,0.01)
nCCs=np.interp(rhos, fiCin, nCCin)
print np.vstack([rhos,nCCs]).T

lCut='3.0'

dinc=1

xn=np.arange(1,2202)-0.5
xp=np.arange(1,2201)

taus=np.zeros([len(rhos)])
nCCsO=np.zeros([len(rhos)])
nCCsO2=np.zeros([len(rhos)])

nMaxS=6
nMaxSMaster=15

for irho in range(len(rhos)):

   fig=plt.figure(num=None, figsize=(6, 6), dpi=600, facecolor='w', edgecolor='w')
   ax = fig.add_subplot(111)

   rho=rhos[irho]
   #print rho
   nCs = [max(20,int(nCCs[irho]-10))]
   nCs = np.append(nCs,np.arange(int(nCCs[irho]-3),int(nCCs[irho]+4)))
   nCs = np.append(nCs,[int(nCCs[irho]+10)])
#   nCs = np.append(nCs,[70])
   taus2=np.zeros([len(nCs)])

   hists = np.zeros([len(xp),len(nCs)])
   ic=0
   for i0 in range(len(nCs)/dinc):
      i=i0*dinc
      nC=nCs[i]
      dirN='nC_'+str(nC)+'_rho_'+'{:.2f}'.format(rho)+'_rC_'+lCut

      #print rho,nC

      counts=np.array([])

      for iF in range(1,4):

         fIn=datRoot+dirN+'/cluster.'+str(iF)+'.dat'
         #print fIn

         res0 = np.fromfile(fIn,sep=" ")
         nRes = len(res0)/nC
         res  = res0.reshape(nC,nRes)
         for iT in range(nC):
            u, count = np.unique(res[iT,:], return_counts=True)
            counts=np.append(counts,count)
      #print counts
      hist=np.histogram(counts, bins=xn, normed=True)

      #np.savetxt('p(s).phi='+str(rho)+'.nC='+str(nC)+'.dat',np.vstack([xp,hist[0]]).T)
      if hist[0][1]>0:
         ax.plot(xp,hist[0],lSs[ic%len(lSs)]+'-' ,label=r'$c=$'+str(nC))
         ic+=1
      hists[:,i] = hist[0]
      np.savetxt('p(s).phi.'+str(rho)+'.nC.'+str(nC)+'.dat',np.vstack([xp,hist[0]]).T)
      
      idxes = hist[0][:nMaxS]>0
      if len(xp[:nMaxS][idxes])>0:
         x = np.log10(xp[:nMaxS][idxes])
         y = np.log10(hist[0][:nMaxS][idxes])
         A = np.vstack([x, np.ones(len(x))]).T
         #print idxes,x,y
         m, c = np.linalg.lstsq(A, y)[0]
         taus2[i0] = np.abs(m)

   taus2[taus2==0]=100
   taus[irho]=np.min(taus2)
   nCCsO[irho]=nCs[ np.argmin(taus2) ]

   c = 1.2
   x = xp[:50]
   y = c*np.power(x, -taus[irho])
   np.savetxt('p(s)asymptotic.phi.'+str(rho)+'.dat',np.vstack([x,y]).T)
   
   ax.plot(x,y,'k--') # ,label=r'\tau= '+'{:}'.format(taus[irho]))
   ax.text(x[3],y[3],r'$\tau = $'+'{:.2f}'.format(taus[irho]),fontsize=20)

   ax.set(xlabel=r'$s$', ylabel=r'$p_{c_{\rm{}}}(s)$',title=r'$\phi=$'+str(rho))
   ax.set_title(r'$\phi=$'+str(rho),fontsize=20)
   #fig.suptitle(r'$\phi=$'+str(rho), fontsize=20)

   legend = ax.legend(loc='upper right', shadow=False, frameon=False, fontsize=13)
   #legend.get_frame().set_facecolor('#00FFCC')
   ax.xaxis.label.set_size(24)
   ax.yaxis.label.set_size(24)
   ax.tick_params(axis='x', labelsize=20)
   ax.tick_params(axis='y', labelsize=20)
   ax.set_xlim(0.9, 50)
   ax.set_ylim(0.000005, 2.)

   for axis in ['top','bottom','left','right']:
       ax.spines[axis].set_linewidth(3.0)

   ax.set_yscale('log')
   ax.set_xscale('log')

   plt.tight_layout()
   
   #legend.get_frame().set_facecolor('#00FFCC')
   figName='clusterSize.phi={:.2f}.png'.format(rho)
   fig.savefig(figName,dpi=600,bbox_inches='tight')
   fig.savefig('clusterSize.phi='+str(rho)+'.svg',dpi=600, format='svg')
   #fig.savefig('clusterSize.phi='+str(rho)+'.eps',dpi=600, format='eps')
   #plt.show()
   plt.close(fig)

   sigmas=np.arange(2.0,2.01,0.1)
   nCMs=np.arange(nCCsO[irho],nCCsO[irho]+0.5,1)
   sqrDists=np.zeros([len(sigmas),len(nCMs)])

   xpM=np.arange(1,13+0.1,0.5)
   histsM  = np.zeros([len(xpM),len(nCs)])
   histsp = np.zeros([len(xp),len(nCs)])
   xsp    = np.zeros([len(xp),len(nCs)])

   isigma=-1
   for sigma in sigmas:
      isigma+=1
      inCM=-1
      for nCM in nCMs:
         inCM+=1
         for i1 in range(len(nCs)):
            #xsp[:,i1]    = xp*np.power(rho*np.abs(27/nCM-27/nCs[i1]),1/sigma)
            xsp[:,i1]    = xp*np.power(np.abs(nCM-nCs[i1]),1/sigma)
            histsp[:,i1] = hists[:,i1]*np.power(xp,taus[irho])
	    histsM[:,i1] = np.interp(xpM, xsp[:,i1], histsp[:,i1])
	    #print np.power(np.abs(27/nCM-27/nCs[i1]),1/sigma),histsM[:,i1]

         for i1 in range(len(nCs)-1):
            for i2 in range(i1+1,len(nCs)):
               distance=np.linalg.norm(np.log10(histsM[:,i1])-np.log10(histsM[:,i2]))
               #print nCs[i1],nCs[i2],distance,sqrDists[isigma,inCM],sigma,nCM
               if distance<1.0 and distance>0:
                  sqrDists[isigma,inCM]+=(distance**2)

   #print sqrDists
   fig=plt.figure(num=None, figsize=(5, 5), dpi=600, facecolor='w', edgecolor='w')
   ax1 = fig.add_subplot(111)
   #ax1 = inset_axes(ax,width="30%",height="30%",loc=3)
   #ax1 = plt.axes([0.3,0.25,0.2,0.2])
   sqrDists[sqrDists==0]=1000000000
   ind = np.unravel_index(np.argmin(sqrDists, axis=None), sqrDists.shape)
   sigma = sigmas[ind[0]]
   nCM   = nCMs[ind[1]]

   nCCsO2[irho]=nCM
   print rho,sigma,nCM,nCCsO[irho]
   ic=0
   for i1 in range(len(nCs)):
      #xsp[:,i1]    = xp*np.power(rho*np.abs(27/nCM-27/nCs[i1]),1/sigma)
      xsp[:,i1]    = xp*np.power(np.abs(nCM-nCs[i1]),1/sigma)
      histsp[:,i1] = hists[:,i1]*np.power(xp,taus[irho])
      if hists[1,i1]>0:
         ax1.plot(xsp[:,i1],histsp[:,i1],lSs[ic%len(lSs)]+'-' ,label=str(nCs[i1]))
         ic+=1

   ax1.set(xlabel=r'$s|c-c^*|^{1/\sigma}$', ylabel=r'$p_{c_{\rm{}}}(s)s^\tau$',title=r'$\phi=$'+str(rho))
   ax1.set_title(r'$\phi=$'+str(rho),fontsize=20)
   #legend = ax.legend(loc='lower left', shadow=False, frameon=False, fontsize=11)
   ax1.xaxis.label.set_size(24)
   ax1.yaxis.label.set_size(24)
   ax1.tick_params(axis='x', labelsize=20)
   ax1.tick_params(axis='y', labelsize=20)
   ax1.set_xlim(0.5, 50)
   ax1.set_ylim(0.01, 2.)

   for axis in ['top','bottom','left','right']:
      ax1.spines[axis].set_linewidth(3.0)

   ax1.set_yscale('log')
   ax1.set_xscale('log')

   #plt.tight_layout()

   figName='clusterSize.master.phi={:.2f}'.format(rho)
   fig.savefig(figName+'.png',dpi=600,bbox_inches='tight')
   fig.savefig(figName+'.svg',dpi=600,format='svg',bbox_inches='tight')
   plt.close(fig)
