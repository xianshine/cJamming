#!/bin/python 

from dump import dump
from pdbfile import pdbfile
import re
import sys
import os
import numpy as np
import scipy.sparse as sp
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import math as m

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

def distance(x0, x1, boxL):
   if len(boxL)==1:
      dimensions=np.ones(len(x0))*boxL
   else :
      dimensions=boxL

   delta = np.abs(x0 - x1)
   delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
   return np.sqrt((delta ** 2).sum(axis=-1))

def pairPot(r,sigma,eps):
   if r<=sigma:
      energy=0.5*eps*(1-r/sigma)**2
   else:
      energy=0
   return energy

dirNow=os.getcwd().split("/")[-1].split("_")
#nB_4200_nC_1_rho_0.60_rC_2.0
#nB=int(dirNow[1])
nC=int(dirNow[1])
rho=float(dirNow[3])
rC=float(dirNow[5])
lcut=rC-5e-5

fcoor='xyz'
fdist='matdist'

#matDist   =sp.csr_matrix((nAperC, nAperC), dtype=np.float32)
maxNumInCluster=np.array([])
numClusters=np.array([])
ratioInMaxCluster=np.array([])

nFrames=3
for i in range(1,nFrames+1):
   dumpName = 'min.'+str(i)+'.atom'
   d = dump(dumpName);
   #p = pdbfile(d);
   #p.one('test.pdb')
   #d.tselect.one("$t >= 100")
   #t = d.time()
   snap=d.snaps[d.nsnaps-1]
   xlo=snap.xlo
   xhi=snap.xhi
   natoms=snap.natoms
   atoms=snap.atoms
   nB=natoms
   nAperC=int(nB/nC)
   box=(xhi-xlo)*np.ones(3)
   res = np.zeros((nC,nAperC),dtype=np.int)-2
   #matDistDense=np.zeros((nAperC, nAperC))

   for iC in range(1,nC+1):
      aiC=atoms[atoms[:,1]==iC]
      np.savetxt(fcoor,aiC[:,2:5],fmt='%20.10e')
      
      #cmdTxt='../../../fileHS2/dist '+str(nAperC) \
      #+' '+str(xlo) \
      #+' '+str(xhi) \
      #+' '+fcoor \
      #+' '+fdist 
      cmdTxt='../../../fileHS2/dist {:d} {:15.10f} {:15.10f} {} {}' \
      .format(nAperC,xlo,xhi,fcoor,fdist)
      #print cmdTxt
      os.system(cmdTxt)
      
      matDistDense=np.fromfile(fdist,sep=" ").reshape(nAperC, nAperC)
      #matDist=sp.csr_matrix(matDistDense)
      #print matDist.todense()
      
      db = DBSCAN(eps=lcut, metric='precomputed',min_samples=1).fit(matDistDense)
      #db = DBSCAN(eps=3.0, metric='precomputed').fit(matDistDense)
      
      labels = db.labels_
      #print db
      #print labels
      res[iC-1,:]=labels
      unique, counts = np.unique(labels, return_counts=True)
      sortCounts=np.sort(counts)
      maxNumInCluster  =np.append(maxNumInCluster,sortCounts[-1])
      numClusters      =np.append(numClusters,len(unique))
      ratioInMaxCluster=np.append(ratioInMaxCluster,float(sortCounts[-1])/float(nAperC))

   f='cluster.'+str(i)+'.dat'
   np.savetxt(f,res,fmt='%4d')
   
sum1=[np.mean(maxNumInCluster.astype(float)),np.std(maxNumInCluster.astype(float))]
sum2=[np.mean(numClusters.astype(float)),np.std(numClusters.astype(float))]
sum3=[np.mean(ratioInMaxCluster.astype(float)),np.std(ratioInMaxCluster.astype(float))]
summary = np.array([sum1,sum2,sum3])
np.savetxt('cluster.summary.dat',summary,fmt='%6.3f')
