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
                    sisj14[i,j,id1,id2] = indat[icount] ##+ (-indat2[icount] if ldeductbg else 0.j)
                    sisj14[j,i,id2,id1] = indat[icount] ##+ (-indat2[icount] if ldeductbg else 0.j)
                    icount += 1
                    sisj23[i,j,id1,id2] = indat[icount] ##+ ( indat2[icount] if ldeductbg else 0.j) #sign canceled
                    sisj23[j,i,id2,id1] = indat[icount] ##+ ( indat2[icount] if ldeductbg else 0.j) #sign canceled
                    icount += 1
            #np.set_printoptions(precision=2,linewidth=400)
            #print( sisj14[i,j,id1,:] )

# calculate ddab (a,b = 1, 2, 3) vs x. Here 1 denotes x, 2 denotes y, 3 denotes xy
# pick a reference point x0 to be x0=Nx/4-1
# id = 0 : +x direc
#      1 : -y direc, fourth quardrant (south-east), -a2 direc.
#      2 : -y direc, third quardrant (south-west)
#      3 : -x direc
#      4 : +y direc, second quardrant (north-west)
#      5 : +y direc, second quardrant (north-east)
dd11_collect = np.zeros((Nx/4,Nx/2),dtype=complex)
dd12_collect = np.zeros((Nx/4,Nx/2),dtype=complex)
dd13_collect = np.zeros((Nx/4,Nx/2),dtype=complex)
dd21_collect = np.zeros((Nx/4,Nx/2),dtype=complex)
dd22_collect = np.zeros((Nx/4,Nx/2),dtype=complex)
dd23_collect = np.zeros((Nx/4,Nx/2),dtype=complex)
dd31_collect = np.zeros((Nx/4,Nx/2),dtype=complex)
dd32_collect = np.zeros((Nx/4,Nx/2),dtype=complex)
dd33_collect = np.zeros((Nx/4,Nx/2),dtype=complex)
for x0 in range(Nx/4-1, Nx/2-1) :
    dd11 = np.zeros((Nx),dtype=complex)
    dd12 = np.zeros((Nx),dtype=complex)
    dd13 = np.zeros((Nx),dtype=complex)
    dd21 = np.zeros((Nx),dtype=complex)
    dd22 = np.zeros((Nx),dtype=complex)
    dd23 = np.zeros((Nx),dtype=complex)
    dd31 = np.zeros((Nx),dtype=complex)
    dd32 = np.zeros((Nx),dtype=complex)
    dd33 = np.zeros((Nx),dtype=complex)
    for j in range(N):
        #if j%Ny==1 :
        i0 = x0*Ny + j%Ny
        dd11[j/Ny] += ( sisj14[i0,j,0,0] - sisj23[i0,j,0,0] )
        dd12[j/Ny] += ( sisj14[i0,j,0,4] - sisj23[i0,j,0,4] )
        dd13[j/Ny] += ( sisj14[i0,j,0,5] - sisj23[i0,j,0,5] )
        dd21[j/Ny] += ( sisj14[i0,j,4,0] - sisj23[i0,j,4,0] )
        dd22[j/Ny] += ( sisj14[i0,j,4,4] - sisj23[i0,j,4,4] )
        dd23[j/Ny] += ( sisj14[i0,j,4,5] - sisj23[i0,j,4,5] )
        dd31[j/Ny] += ( sisj14[i0,j,5,0] - sisj23[i0,j,5,0] )
        dd32[j/Ny] += ( sisj14[i0,j,5,4] - sisj23[i0,j,5,4] )
        dd33[j/Ny] += ( sisj14[i0,j,5,5] - sisj23[i0,j,5,5] )
    dd11_collect[x0-(Nx/4-1)] = dd11[x0:x0+Nx/2]
    dd12_collect[x0-(Nx/4-1)] = dd12[x0:x0+Nx/2]
    dd13_collect[x0-(Nx/4-1)] = dd13[x0:x0+Nx/2]
    dd21_collect[x0-(Nx/4-1)] = dd21[x0:x0+Nx/2]
    dd22_collect[x0-(Nx/4-1)] = dd22[x0:x0+Nx/2]
    dd23_collect[x0-(Nx/4-1)] = dd23[x0:x0+Nx/2]
    dd31_collect[x0-(Nx/4-1)] = dd31[x0:x0+Nx/2]
    dd32_collect[x0-(Nx/4-1)] = dd32[x0:x0+Nx/2]
    dd33_collect[x0-(Nx/4-1)] = dd33[x0:x0+Nx/2]
    
    with open("dd11vsx_x0"+str(x0)+".dat","w") as f:
        for ix in range(x0,Nx):
            f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd11[ix].real, dd11[ix].imag) )
    with open("dd12vsx_x0"+str(x0)+".dat","w") as f:
        for ix in range(x0,Nx):
            f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd12[ix].real, dd12[ix].imag) )
    with open("dd13vsx_x0"+str(x0)+".dat","w") as f:
        for ix in range(x0,Nx):
            f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd13[ix].real, dd13[ix].imag) )
    with open("dd21vsx_x0"+str(x0)+".dat","w") as f:
        for ix in range(x0,Nx):
            f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd21[ix].real, dd21[ix].imag) )
    with open("dd22vsx_x0"+str(x0)+".dat","w") as f:
        for ix in range(x0,Nx):
            f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd22[ix].real, dd22[ix].imag) )
    with open("dd23vsx_x0"+str(x0)+".dat","w") as f:
        for ix in range(x0,Nx):
            f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd23[ix].real, dd23[ix].imag) )
    with open("dd31vsx_x0"+str(x0)+".dat","w") as f:
        for ix in range(x0,Nx):
            f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd31[ix].real, dd31[ix].imag) )
    with open("dd32vsx_x0"+str(x0)+".dat","w") as f:
        for ix in range(x0,Nx):
            f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd32[ix].real, dd32[ix].imag) )
    with open("dd33vsx_x0"+str(x0)+".dat","w") as f:
        for ix in range(x0,Nx):
            f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd33[ix].real, dd33[ix].imag) )
dd11 = np.mean(dd11_collect, axis=0)
dd12 = np.mean(dd12_collect, axis=0)
dd13 = np.mean(dd13_collect, axis=0)
dd21 = np.mean(dd21_collect, axis=0)
dd22 = np.mean(dd22_collect, axis=0)
dd23 = np.mean(dd23_collect, axis=0)
dd31 = np.mean(dd31_collect, axis=0)
dd32 = np.mean(dd32_collect, axis=0)
dd33 = np.mean(dd33_collect, axis=0)

with open("dd11vsx.dat","w") as f:
    for ix in range(len(dd11)):
        f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd11[ix].real, dd11[ix].imag) )
with open("dd12vsx.dat","w") as f:
    for ix in range(len(dd11)):
        f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd12[ix].real, dd12[ix].imag) )
with open("dd13vsx.dat","w") as f:
    for ix in range(len(dd11)):
        f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd13[ix].real, dd13[ix].imag) )
with open("dd21vsx.dat","w") as f:
    for ix in range(len(dd11)):
        f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd21[ix].real, dd21[ix].imag) )
with open("dd22vsx.dat","w") as f:
    for ix in range(len(dd11)):
        f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd22[ix].real, dd22[ix].imag) )
with open("dd23vsx.dat","w") as f:
    for ix in range(len(dd11)):
        f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd23[ix].real, dd23[ix].imag) )
with open("dd31vsx.dat","w") as f:
    for ix in range(len(dd11)):
        f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd31[ix].real, dd31[ix].imag) )
with open("dd32vsx.dat","w") as f:
    for ix in range(len(dd11)):
        f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd32[ix].real, dd32[ix].imag) )
with open("dd33vsx.dat","w") as f:
    for ix in range(len(dd11)):
        f.write( "{:4d} {: .12f} {: .12f}\n".format(ix, dd33[ix].real, dd33[ix].imag) )
