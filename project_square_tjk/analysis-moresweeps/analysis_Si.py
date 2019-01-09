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
assert len(sys.argv) == 3, "Usage: python file.py datafile tag"
tag=sys.argv[2]+"sp"

# read data
print "reading file "+sys.argv[1]+" ......"
indati=np.loadtxt(sys.argv[1])  #
numsi = len(indati)
print "len(indati) = ", numsi
bg=np.mean(indati)

## check Si.out, to take care different cases
yfold = 1 # a factor for taking care of chiral correlation
if yperiodic :
    if numsi == N:
        Neff = N  ## for spin correlation, y_dimer correlation
    elif numsi == (N-Ny):
        Neff = N-Ny ## for x_dimer, xy_dimer correlation
    elif numsi == (N-Ny)*2:
        Neff = 2*(N-Ny) ## for chiral correlation
        yfold = 2 ## two tri_plaq for one unit cell, we count it in y direction
    else:
        print "Wrong number of SiSj, please check your SiSj.out"
elif Ny > 2: # y-direction open boundary case, and not the ladder
    if numsi == N:
        Neff = N  ## for spin correlation
    elif numsi == (N-Ny):
        Neff = N-Ny ## for x_dimer correlation
    elif numsi == (N-Nx):
        Neff = N-Nx ## for y_dimer correlation
    elif numsi == (N-Nx-Ny+1):
        Neff = N-Nx-Ny+1 ## for xy_dimer correlation
    elif numsi == (N-Nx-Ny+1)*2:
        Neff = 2*(N-Nx-Ny+1) ## for chiral correlation
        yfold = 2 ## two tri_plaq for one unit cell, we count it in y direction
    else:
        print "Wrong number of SiSj, please check your SiSj.out"
else : # ladder case
    if numsi == N:
        Neff = N  ## for spin correlation
    elif numsi == (N-Ny):
        Neff = N-Ny ## for x_dimer correlation
    elif numsi == (N-Nx):
        Neff = N-Nx ## for y_dimer correlation
    elif numsi == (N-Nx-Ny+1):
        Neff = N-Nx-Ny+1 ## for xy_dimer correlation
    elif numsi == (N-Nx-Ny+1)*2:
        Neff = 2*(N-Nx-Ny+1) ## for chiral correlation
    else:
        print "Wrong number of SiSj, please check your SiSj.out"

kfac = 1
## redefine N and Ny
N *= yfold
Ny *= yfold
if Ny == 2 and (not yperiodic) :
    Ny /= 2
    Nx *= 2
# k point coordinate, consider full boundary for better plot
    Nk = (Nx+1)*Ny
    kxv = np.zeros(Nk)
    kyv = np.zeros(Nk)
    for j in range(Ny):
        for i in range(Nx+1):
            kxv[(Nx+1)*j+i] = float(i) - float(Nx)/2.0
            kyv[(Nx+1)*j+i] = 0.0
else :
# k point coordinate, consider full boundary for better plot
    Nk = (kfac*Nx+1)*(kfac*Ny+1)
    kxv = np.zeros(Nk)
    kyv = np.zeros(Nk)
    for j in range(kfac*Ny+1):
        for i in range(kfac*Nx+1):
            kxv[(kfac*Nx+1)*j+i] = float(i) - float(kfac*Nx)/2.0
            kyv[(kfac*Nx+1)*j+i] = float(j) - float(kfac*Ny)/2.0

# real space coordinate
xv = np.zeros(N)
yv = np.zeros(N)
for i in range(Nx):
    for j in range(Ny):
        xv[Ny*i+j] = float(i)
        yv[Ny*i+j] = float(j)

# set fourier phase
expirk = np.zeros((N,Nk),dtype=complex)
for ri in range(N):
    for kj in range(Nk):
        expirk[ri,kj] = np.exp( 1.j*(xv[ri]*kxv[kj]/Nx+yv[ri]*kyv[kj]/Ny)*2.0*np.pi )
        ##print "expirk[ri,kj]=", expirk[ri,kj]

simj = np.zeros(N)
# load indati into simj, deduct background
for i in range(2*Ny,len(indati)-2*Ny):  # Note boundary data are dropped out
    simj[i] = indati[i]-bg

# perform fourier transformation
with open(tag+"k.dat","w") as f:
    for ki in range(Nk):
        ssk = 0.+0.j
        for rimj in range(N):
            ssk += simj[rimj]*expirk[rimj,ki]
        ssk /= N
        if ki%(kfac*Nx+1) == 0 and ki!=0 :
            f.write( "\n" )
        #f.write( "{: .8f} {: .8f} {: .8f} {: .8f}\n".format(kxv[ki]/Nx*2.0*np.pi, kyv[ki]/Ny*2.0*np.pi, ssk.real, ssk.imag) )
        f.write( "{: .8f} {: .8f} {: .8f} {: .8f}\n".format(kxv[ki]/Nx, kyv[ki]/Ny, ssk.real, ssk.imag) )
