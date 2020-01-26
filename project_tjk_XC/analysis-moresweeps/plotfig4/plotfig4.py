#!/usr/bin/dev python
#coding=utf-8
"""
Author:         Xiao-Yan Xu <wanderxu@gmail.com>
Description:
use chi square to do fitting.

Input:
      file "tmppy.dat" with first k columns of independent variables,
                            following columns with observables and errors
      take care of the number of independent variables k and the model to fitting defined in curvefunc
      for any special case, you need change corresponding places

"""
import sys
import math
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import ScalarFormatter

## define read file func
def file2list(filename):
		fr = open(filename)
		array = fr.readlines()
		num = len(array)
		returnMat=np.zeros((num,4))  # you can change the dimension 
		index = 0
		for line in array:
				line = line.strip()
				linelist = line.split()
				returnMat[index,:]=linelist[0:4] # you can change the dimension
				index +=1
		return returnMat

## define curvefunc for curve_fit
def curvefunc (xv, *p0 ):
        results = p0[0]+np.log(xv)*p0[1]
        return results

def chi_square ( xdata, ydata, ydata_sigma, p0 ):
        popt,pcov = opt.curve_fit(curvefunc, xdata, ydata, p0, sigma=ydata_sigma, absolute_sigma=False )
        perr = np.sqrt(np.diag(pcov))
        rchi_sq = np.sum( ( (ydata-curvefunc(xdata, *popt ) ) / ydata_sigma )**2 ) / len(ydata)
        return popt, perr, rchi_sq

#assert len(sys.argv) == 2, "Usage: python file.py data.out"
#print "reading file "+sys.argv[1]+" ......"
#indat=np.transpose(np.loadtxt(sys.argv[1]))  #

nfile=6
x1=[None]*nfile
y1=[None]*nfile
dy=[None]*nfile
for ifile in range(nfile):
    indat=np.transpose(np.loadtxt("file"+"%d"%ifile+".dat"))
    ####orignial data
    x1[ifile]=indat[0]
    y1[ifile]=indat[2]
    dy[ifile] = np.zeros(y1[ifile].shape,np.float)
    dy[ifile][0:-1] = np.diff(y1[ifile])/np.diff(x1[ifile])
    dy[ifile][-1] = (y1[ifile][-1] - y1[ifile][-2])/(x1[ifile][-1] - x1[ifile][-2])
    
fig = plt.figure(figsize=(8, 6)) 
nrow=3
ncol=2
ax=[[None]*ncol]*nrow
#print ax[0]
kytag=['$k_2=-0.25$\n$\\frac{1}{8}$ dop',
       '$k_2=-0.25$\n$\\frac{1}{12}$ dop',
       '$k_2=-0.25$\n$\\frac{1}{16}$ dop',
       '$k_2=0$\n$\\frac{1}{8}$ dop',
       '$k_2=0$\n$\\frac{1}{12}$ dop',
       '$k_2=0$\n$\\frac{1}{16}$ dop']
figtag=['a','c','e','b','d','f']
#leftp=[-3.0/8, -5.0/12, -5.0/12, -7.0/24, -7.0/24, -7.0/24]
#rightp=[1.0/4, 1.0/4, 1.0/4, 7.0/24, 7.0/24, 7.0/24]
leftp=[-1.0/4, -1.0/4, -1.0/4, -7.0/24, -7.0/24, -7.0/24]
rightp=[3.0/8, 5.0/12, 5.0/12, 7.0/24, 7.0/24, 7.0/24]
gs = gridspec.GridSpec(nrow, ncol, height_ratios=[1, 1, 1]) 
ifile = 0
for icol in range(ncol):
    for irow in range(nrow):
        ax[irow][icol] = plt.subplot(gs[irow,icol])
        ax[irow][icol].plot(x1[ifile], y1[ifile], '.b', label=kytag[ifile])
        ax[irow][icol].set_ylabel('$n(\mathbf{k})$',color='b')
        ax[irow][icol].tick_params('y', colors='b')

        major_ticks = np.arange(-0.5, 0.51, 0.5)
        minor_ticks = np.arange(-0.5, 0.51, 1.0/12)
        ax[irow][icol].set_xticks(major_ticks)
        ax[irow][icol].set_xticks(minor_ticks, minor=True)
        ax[irow][icol].grid(which='both', axis='x')
        ax[irow][icol].text(-0.16, 1.05, figtag[ifile], transform=ax[irow][icol].transAxes, fontsize=14, fontweight='bold', va='top', ha='right')
        ax[irow][icol].legend(fontsize=10)

        axd = ax[irow][icol].twinx()
        axd.plot(x1[ifile], dy[ifile], '-r')
        axd.set_ylabel('$\partial n(\mathbf{k}) / \partial k_1 $',color='r')
        axd.tick_params('y', colors='r')
        #axd.plot([leftp[ifile],leftp[ifile]], [min(dy[ifile]),max(dy[ifile])], '--r') # the peak position
        #axd.plot([rightp[ifile],rightp[ifile]], [min(dy[ifile]),max(dy[ifile])], '--r') # the peak position
        dr=(max(dy[ifile])-min(dy[ifile]))*0.05
        axd.arrow(leftp[ifile],  max(dy[ifile])-6*dr, 0.0,  2*dr, color='k', head_length=1*dr, head_width=0.02, linewidth=0.1 ) # the peak position
        axd.arrow(rightp[ifile], min(dy[ifile])+6*dr, 0.0, -2*dr, color='k', head_length=1*dr, head_width=0.02, linewidth=0.1 ) # the peak position
        ifile = ifile + 1

plt.tight_layout()
plt.savefig("fig4.pdf", dpi=300)
