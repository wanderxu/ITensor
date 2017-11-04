#!/usr/bin/dev python
#coding=utf-8
#
# can be used to analysis spin, dimer and chiral correlation
#

import sys
import math
import numpy as np
import scipy.optimize as opt

# set parameter
from model_para import *

assert len(sys.argv) == 3, "Usage: python file.py data.out tag"
tag=sys.argv[2]

# read data
print "reading file "+sys.argv[1]+" ......"
indat=np.loadtxt(sys.argv[1])  #
numsij = len(indat)
print "len(indat) = ", numsij

## check SiSj.out, to take care different cases
yfold = 1 # a factor for taking care of chiral correlation
if numsij == N*(N+1)/2:
    Neff = N  ## for spin correaltion, y_dimer correaltion
elif numsij == (N-Ny)*(N-Ny+1)/2:
    Neff = N-Ny ## for x_dimer, xy_dimer correaltion
elif numsij == (N-Ny)*(2*N-2*Ny+1):
    Neff = 2*(N-Ny) ## for chiral correlation
    yfold = 2 ## two tri_plaq for one unit cell, we count it in y direction
else:
    print "Wrong number of SiSj, please check your SiSj.out"

## redefine N and Ny
N *= yfold
Ny *= yfold

# set fourier phase
expirk = np.zeros((N,N),dtype=complex)
xi = np.zeros(N)
yi = np.zeros(N)
for i in range(Nx):
    for j in range(Ny):
        xi[Ny*i+j] = float(i)
        yi[Ny*i+j] = float(j)
kxi = xi - float(Nx)/2.0
kyi = yi - float(Ny)/2.0

pra = np.array([1.0,0.0])
prb = np.array([-0.5,np.sqrt(3.0)/2.0])
pka = np.array([np.sqrt(3.0)/2.0,0.5])*4.0*np.pi/np.sqrt(3.0)/float(Nx)
pkb = np.array([0.0,-1.0])            *4.0*np.pi/np.sqrt(3.0)/float(Ny)

kk = np.zeros((N,2))
rr = np.zeros((N,2))
for i in range(N):
    kk[i][:] = kxi[i]*pka + kyi[i]*pkb
    rr[i][:] = xi[i]*pra + yi[i]*prb

for ri in range(N):
    for kj in range(N):
        expirk[ri,kj] = np.exp( 1.j*np.inner(rr[ri],kk[kj]) )
        #print "expirk[ri",ri,"kj",kj,"]=", expirk[ri,kj]

sisj = np.zeros((N,N))
simj = np.zeros(N)

# load indat to sisj
icount = 0
for i in range(Neff):
    for j in range(i,Neff):
        if icount < numsij:
            sisj[i,j] = indat[icount]
            sisj[j,i] = indat[icount]
            icount += 1
    np.set_printoptions(precision=2,linewidth=400)
    print( sisj[i][0:max(Neff,20)] )
    #print str(sisj[i])
    #print( sisj[i], end=" ")
    #print sisj[i],
    #print "16{: .8f}".format(sisj[i])

# calculate simj
for i in range(N):
    for j in range(N):
        imj = int( (xv[i]-xv[j] + Nx)%Nx * Ny + (yv[i]-yv[j]+Ny)%Ny )
        #print "imj=",imj
        simj[imj] += sisj[i,j]

# perform fourier transformation
with open(tag+"k.dat","w") as f:
    for ki in range(N):
        ssk = 0.+0.j
        for rimj in range(N):
            ssk += simj[rimj]*expirk[rimj,ki]
        ssk /= float(N)
        if ki%Ny == 0 :
            f.write( "\n" )
        f.write( "{: .8f} {: .8f} {: .8f} {: .8f}\n".format(kk[ki][0], kk[ki][1], ssk.real, ssk.imag) )

# output sisj in x-direction
# pick the central point
ic=(Nx/2-1)*Ny+Ny/2-1

# x-direction
with open(tag+"ij_xdirec.dat","w") as f:
    for i in range(ic,N,Ny):
        f.write( "{} {: .8f}\n".format((i/Ny-Nx/2+1)*yfold, sisj[ic][i]) )
        if yfold == 2 :
            ## in chiral case, you can also count plaq in x direction 
            f.write( "{} {: .8f}\n".format((i/Ny-Nx/2+1)*yfold+1, sisj[ic][i-1]) )

# y-direction
with open(tag+"ij_ydirec.dat","w") as f:
    for i in range(ic,ic+Ny/2+1):
        f.write( "{} {: .8f}\n".format(i-ic, sisj[ic][i]) )
