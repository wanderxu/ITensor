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

# data1.out is correlation data, data2.out is background
assert len(sys.argv) == 4, "Usage: python file.py data1.out data2.out tag"
tag=sys.argv[3]+"dbg"

# read data
print "reading file "+sys.argv[1]+" ......"
indat=np.loadtxt(sys.argv[1])  #
numsij = len(indat)
print "len(indat) = ", numsij

print "reading file "+sys.argv[2]+" ......"
indati=np.loadtxt(sys.argv[2])  #
numsi = len(indati)
print "len(indati) = ", numsi

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
xv = np.zeros(N)
yv = np.zeros(N)
for i in range(Nx):
    for j in range(Ny):
        xv[Ny*i+j] = float(i)
        yv[Ny*i+j] = float(j)
kxv = xv - Nx/2
kyv = yv - Ny/2

for ri in range(N):
    for kj in range(N):
        expirk[ri,kj] = np.exp( 1.j*(xv[ri]*kxv[kj]/Nx+yv[ri]*kyv[kj]/Ny)*2.0*np.pi )
        ##print "expirk[ri,kj]=", expirk[ri,kj]

sisj = np.zeros((N,N))
simj = np.zeros(N)

# load indat to sisj
icount = 0
for i in range(Neff):
    for j in range(i,Neff):
        if icount < numsij:
            sisj[i,j] = indat[icount] - indati[i]*indati[j]
            sisj[j,i] = indat[icount] - indati[j]*indati[i]
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
        ssk /= N
        if ki%Ny == 0 :
            f.write( "\n" )
        #f.write( "{: .8f} {: .8f} {: .8f} {: .8f}\n".format(kxv[ki]/Nx*2.0*np.pi, kyv[ki]/Ny*2.0*np.pi, ssk.real, ssk.imag) )
        f.write( "{: .8f} {: .8f} {: .8f} {: .8f}\n".format(kxv[ki]/Nx, kyv[ki]/Ny, ssk.real, ssk.imag) )

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
