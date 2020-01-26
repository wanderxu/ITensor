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

# data.out is correlation data
assert (len(sys.argv) == 3 or len(sys.argv) == 4), "Usage: python file.py data.out [data2.out] tag"
ldeductbg=(len(sys.argv) == 4)
tag=sys.argv[2]
if ldeductbg:
    tag=sys.argv[3]

# read data.out
print "reading file "+sys.argv[1]+" ......"
indat=np.loadtxt(sys.argv[1], converters={0: lambda s: eval(s)}).view(complex).reshape(-1)
numsij = len(indat)
print "len(indat) = ", numsij

# read data2.out
indat2=0
if ldeductbg:
    print "reading file "+sys.argv[2]+" ......"
    indat2=np.loadtxt(sys.argv[2], converters={0: lambda s: eval(s)}).view(complex).reshape(-1)
    print "len(indat2) = ", len(indat2)

kfac = 1
## redefine N and Ny
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

sisj14 = np.zeros((N,N,6,6),dtype=complex)
sisj23 = np.zeros((N,N,6,6),dtype=complex)
simj = np.zeros(N)

# load indat to sisj14 and sisj23, and calculate <s^2>
icount = 0
for i in range(N):
    for id1 in range(6):
        for j in range(i,N):
            for id2 in range(6):
                if icount < numsij:
                    sisj14[i,j,id1,id2] = indat[icount] + (-indat2[icount] if ldeductbg else 0.j)
                    sisj14[j,i,id2,id1] = indat[icount] + (-indat2[icount] if ldeductbg else 0.j)
                    icount += 1
                    sisj23[i,j,id1,id2] = indat[icount] + ( indat2[icount] if ldeductbg else 0.j) #sign canceled
                    sisj23[j,i,id2,id1] = indat[icount] + ( indat2[icount] if ldeductbg else 0.j) #sign canceled
                    icount += 1
            #np.set_printoptions(precision=2,linewidth=400)
            #print( sisj14[i,j,id1,:] )

##calculate pairing
phase = np.zeros((14,6),dtype=complex)
phase[0,:] = 6*[1.0+0.j]/np.sqrt(6.0)
phase[1,:] = [1.0+0.j, np.exp(2.0*np.pi/3.0*1.j), np.exp(-2.0*np.pi/3.0*1.j), 1.0+0.j, np.exp(2.0*np.pi/3.0*1.j), np.exp(-2.0*np.pi/3.0*1.j)]/np.sqrt(6.0)
phase[2,:] = [1.0+0.j, np.exp(-2.0*np.pi/3.0*1.j), np.exp(2.0*np.pi/3.0*1.j), 1.0+0.j, np.exp(-2.0*np.pi/3.0*1.j), np.exp(2.0*np.pi/3.0*1.j)]/np.sqrt(6.0)
phase[3,:] = [0.0+0.j, -1.0+0.j, 1.0+0.j, 0.0+0.j, -1.0+0.j, 1.0+0.j]/np.sqrt(4.0)
phase[4,:] = [2.0+0.j, -1.0+0.j, -1.0+0.j, 2.0+0.j, -1.0+0.j, -1.0+0.j]/np.sqrt(12.0)
phase[5,:] = [1.0+0.j, -1.0+0.j, 1.0+0.j, -1.0+0.j, 1.0+0.j, -1.0+0.j]/np.sqrt(6.0)
phase[6,:] = [1.0+0.j, np.exp(np.pi/3.0*1.j), np.exp(2.0*np.pi/3.0*1.j), -1.0+0.j, -np.exp(np.pi/3.0*1.j), -np.exp(2.0*np.pi/3.0*1.j)]/np.sqrt(6.0)
phase[7,:] = [1.0+0.j, np.exp(-np.pi/3.0*1.j), np.exp(-2.0*np.pi/3.0*1.j), -1.0+0.j, -np.exp(-np.pi/3.0*1.j), -np.exp(-2.0*np.pi/3.0*1.j)]/np.sqrt(6.0)
phase[8,:] = [2.0+0.j, 1.0+0.j, -1.0+0.j, -2.0+0.j, -1.0+0.j, 1.0+0.j]/np.sqrt(12.0)
phase[9,:] = [0.0+0.j, -1.0+0.j, -1.0+0.j, 0.0+0.j, 1.0+0.j, 1.0+0.j]/np.sqrt(4.0)
np.set_printoptions(precision=4,linewidth=400)
pairlist=['s', 'd+id', 'd-id', 'dxy', 'dx2-y2', 'f', 'p+ip', 'p-ip', 'px', 'py', '11','12','21','22']
print "phi(", pairlist, ") = "
#print "phi(s, d+id, d-id, dxy, dx2-y2, f, p+ip, p-ip, px, py) = "
print phase

# calculate didj
# drop out boundaries
didj = np.zeros((14,N,N),dtype=complex)
for id1 in range(6):
    for id2 in range(6):
        for i in range(2*Ny,N-2*Ny): # note boundary data are dropped out
            for j in range(2*Ny,N-2*Ny):
                # singlet pairing
                for ipair in range(5):
                    didj[ipair,i,j] += np.conj(phase[ipair,id1])*phase[ipair,id2]*(sisj14[i,j,id1,id2] - sisj23[i,j,id1,id2])
                for ipair in range(5,10):
                # triplet pairing
                    didj[ipair,i,j] += np.conj(phase[ipair,id1])*phase[ipair,id2]*(sisj14[i,j,id1,id2] + sisj23[i,j,id1,id2])
for i in range(2*Ny, N-2*Ny):
    for j in range(2*Ny,N-2*Ny):
        # singlet pairing
        didj[10,i,j] += (sisj14[i,j,0,0] - sisj23[i,j,0,0])
        didj[11,i,j] += (sisj14[i,j,0,4] - sisj23[i,j,0,4])
        didj[12,i,j] += (sisj14[i,j,4,0] - sisj23[i,j,4,0])
        didj[13,i,j] += (sisj14[i,j,4,4] - sisj23[i,j,4,4])

# calculate dimj
# multiply a power to be able to see features in structure factor as it decays so fast
dimj = np.zeros((14,N),dtype=complex)
for i in range(N):
    for j in range(N):
        # as we use PBC in y direction, the length calculation will change a little bit.
        rlen_power = ( (xv[i]-xv[j])**2 + ( float(Ny)-abs(yv[i]-yv[j]) if abs(yv[i]-yv[j]) > float(Ny)/2.0 else abs(yv[i]-yv[j]) )**2 )**3 # r^6
        imj = int( (xv[i]-xv[j] + Nx)%Nx * Ny + (yv[i]-yv[j]+Ny)%Ny )
        #print "imj=",imj
        for ipair in range(14):
            dimj[ipair,imj] += didj[ipair,i,j]*rlen_power

# perform fourier transformation
for ipair in range(14):
    with open(tag+pairlist[ipair]+"k.dat","w") as f:
        for ki in range(Nk):
            ssk = 0.+0.j
            for rimj in range(N):
                ssk += dimj[ipair,rimj]*expirk[rimj,ki]
            ssk /= N
            if ki%(kfac*Nx+1) == 0 and ki!=0 :
                f.write( "\n" )
            #f.write( "{: .8f} {: .8f} {: .8f} {: .8f}\n".format(kxv[ki]/Nx*2.0*np.pi, kyv[ki]/Ny*2.0*np.pi, ssk.real, ssk.imag) )
            f.write( "{: .8f} {: .8f} {: .8f} {: .8f}\n".format(kxv[ki]/Nx, kyv[ki]/Ny, ssk.real, ssk.imag) )

#### output sisj in x-direction
#### pick the central point
###ic=(Nx/2-1)*Ny+Ny/2-1
###
#### x-direction
###with open(tag+"ij_xdirec.dat","w") as f:
###    for i in range(ic,N,Ny):
###        f.write( "{} {: .8f}\n".format((i/Ny-Nx/2+1), sisj[ic][i]) )
###
#### y-direction
###with open(tag+"ij_ydirec.dat","w") as f:
###    for i in range(ic,ic+Ny/2+1):
###        f.write( "{} {: .8f}\n".format(i-ic, sisj[ic][i]) )
