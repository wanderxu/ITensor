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
        #results = p0[0]*xv**p0[1]
        results = p0[0]
        for i in range(len(p0)-1):
            results = results + p0[i+1]*xv**(i+1)
        return results

def chi_square ( xdata, ydata, ydata_sigma, p0 ):
        popt,pcov = opt.curve_fit(curvefunc, xdata, ydata, p0, sigma=ydata_sigma, absolute_sigma=False )
        perr = np.sqrt(np.diag(pcov))
        rchi_sq = np.sum( ( (ydata-curvefunc(xdata, *popt ) ) / ydata_sigma )**2 ) / len(ydata)
        return popt, perr, rchi_sq

assert len(sys.argv) == 2, "Usage: python file.py data.out"
#print "reading file "+sys.argv[1]+" ......"
indat=np.transpose(np.loadtxt(sys.argv[1]))  #

### full points
##xdata = indat[0]
##ydata = indat[1]
##ydata_sigma = len(xdata)*[1]
##p0=4*[1.0] # third order
##popt,perr,rchi_sq = chi_square( xdata, ydata, ydata_sigma, p0 )
##print "full points, 3order, theta0 = ", popt[0], "+/-", perr[0], " x^2 = ", rchi_sq

# drop several points
idrop=3
xdata = indat[0][idrop:]
ydata = np.zeros(len(xdata))
for i in range(len(xdata)):
    if np.abs(indat[1][idrop+i]) > 0:
        ydata[i] = np.log(np.abs(indat[1][idrop+i]))
    else :
        break
ndata = i+1
xdata = xdata[range(ndata)]
ydata = ydata[range(ndata)]
ydata_sigma = len(xdata)*[1]
#p0=2*[1.0] #
#popt,perr,rchi_sq = chi_square( xdata, ydata, ydata_sigma, p0 )
#print "drop 1point, fitting results :"
#print "a = ", popt[0], "+/-", perr[0]
#print "b = ", popt[1], "+/-", perr[1]
#print " x^2 = ", rchi_sq

p0=6*[1.0] #
popt,perr,rchi_sq = chi_square( xdata, ydata, ydata_sigma, p0 )
#print "drop 1point, fitting results :"
for i in range(len(p0)):
    print popt[i], "*x**", i
#print " x^2 = ", rchi_sq

# normalization
ydata_fit = curvefunc(xdata, *popt)
ydata_norm = indat[1][idrop:ndata+idrop]*np.exp(-ydata_fit)

# fft
nk = 200
sp = np.fft.fft(ydata_norm, n=nk) / nk
freq = np.fft.fftfreq(nk)
sp = sp[range(nk/2)]
freq = freq[range(nk/2)]

minor_ticks=np.arange(10)/20.0

fig, ax = plt.subplots(2, 1)
ax[0].plot(xdata,ydata_norm, 'o-')
ax[0].set_xlabel('$x-x_0$', fontsize=12)
ax[0].set_ylabel('Normalized correlation')
ax[0].axhline(y=0, color='k', linestyle='--', linewidth=0.5)
ax[1].plot(freq,abs(sp),'r') # plotting the spectrum
ax[1].set_xlabel('$Q/2\pi$')
ax[1].set_ylabel('|F(Q)|')
ax[1].set_xticks(minor_ticks, minor=True)
ax[1].grid(which='major')
ax[1].grid(which='minor')

plt.tight_layout()
plt.savefig("tmp.pdf", dpi=300)
