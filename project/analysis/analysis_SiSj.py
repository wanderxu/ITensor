#!/usr/bin/dev python
#coding=utf-8

import sys
import math
import numpy as np
import scipy.optimize as opt

# set parameter
from model_para import *

assert len(sys.argv) == 2, "Usage: python file.py data.out"

# read data
print "reading file "+sys.argv[1]+" ......"
indat=np.loadtxt(sys.argv[1])  #
#print "input = ", indat

assert len(indat) == N*(N+1)/2, "Wrong number of SiSj, please check your SiSj.out"

# set fourier phase
expirk = np.zeros((N,N),dtype=complex)
xv = np.zeros(N)
yv = np.zeros(N)
for i in range(Nx):
    for j in range(Ny):
        xv[Ny*i+j] = float(i)
        yv[Ny*i+j] = float(j)

for ri in range(N):
    for kj in range(N):
        expirk[ri,kj] = np.exp( 1.j*(xv[ri]*xv[kj]/Nx+yv[ri]*yv[kj]/Ny)*2.0*np.pi )
        ##print "expirk[ri,kj]=", expirk[ri,kj]

sisj = np.zeros((N,N))
simj = np.zeros(N)

# load indat to sisj and simj
icount = 0
for i in range(N):
    for j in range(i,N):
        sisj[i,j] = indat[icount]
        sisj[j,i] = indat[icount]
        imj = int( (xv[i]-xv[j] + Nx)%Nx * Ny + (yv[i]-yv[j]+Ny)%Ny )
        #print "imj=",imj
        simj[imj] += indat[icount]
        icount += 1

# perform fourier transformation
with open("ssk.dat","w") as f:
    for ki in range(N):
        ssk = 0.+0.j
        for rimj in range(N):
            ssk += simj[rimj]*expirk[rimj,ki]
        ssk /= N
        f.write( "{: .8f} {: .8f} {: .8f} {: .8f}\n".format(xv[ki]/Nx*2.0*np.pi, yv[ki]/Ny*2.0*np.pi, ssk.real, ssk.imag) )

# output sisj in x-direction
# pick the central point
ic=(Nx/2-1)*Ny+Ny/2-1

# x-direction
with open("ssij_xdirec.dat","w") as f:
    for i in range(ic,N,Ny):
        f.write( "{} {: .8f}\n".format(i/Nx-Nx/2+1, sisj[ic][i]) )

# y-direction
with open("ssij_ydirec.dat","w") as f:
    for i in range(ic,ic+Ny/2+1):
        f.write( "{} {: .8f}\n".format(i-ic, sisj[ic][i]) )
