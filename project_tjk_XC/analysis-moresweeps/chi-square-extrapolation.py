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
        results = p0[0]
        for i in range(len(p0)-1):
            results = results + p0[i+1]*xv**(i+1)
        return results

def chi_square ( xdata, ydata, ydata_sigma, p0 ):
        popt,pcov = opt.curve_fit(curvefunc, xdata, ydata, p0, sigma=ydata_sigma, absolute_sigma=False )
        perr = np.sqrt(np.diag(pcov))
        rchi_sq = np.sum( ( (ydata-curvefunc(xdata, *popt ) ) / ydata_sigma )**2 ) / len(ydata)
        return popt, perr, rchi_sq

assert len(sys.argv) == 4, "Usage: python file.py data1.out data2.out data3.out"
#print "reading file "+sys.argv[1]+" ......"
indat1=np.transpose(np.loadtxt(sys.argv[1]))  #
indat2=np.transpose(np.loadtxt(sys.argv[2]))  #
indat3=np.transpose(np.loadtxt(sys.argv[3]))  #

xdata = indat1[0]
ydata1 = indat1[1]
ydata2 = indat2[1]
ydata3 = indat3[1]
xi=1.0/np.array([1280, 2560, 5120])
yi_sigma=3*[1]
p0=3*[1.0] #
for i in range(len(ydata1)):
    yi=np.array([ydata1[i], ydata2[i], ydata3[i]])
    popt,perr,rchi_sq = chi_square( xi, yi, yi_sigma, p0 )
    print xdata[i], popt[0]
