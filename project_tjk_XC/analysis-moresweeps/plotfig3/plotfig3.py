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
xdata=[None]*nfile
ydata_norm=[None]*nfile
popt=[None]*nfile
pos_y=[None]*nfile
neg_y=[None]*nfile
sp=[None]*nfile
freq=[None]*nfile
x1=[None]*nfile
y1=[None]*nfile
for ifile in range(nfile):
    indat=np.transpose(np.loadtxt("file"+"%d"%ifile+".dat"))
    print indat

    ### full points
    ##xdata = indat[0]
    ##ydata = indat[1]
    ##ydata_sigma = len(xdata)*[1]
    ##p0=4*[1.0] # third order
    ##popt,perr,rchi_sq = chi_square( xdata, ydata, ydata_sigma, p0 )
    ##print "full points, 3order, theta0 = ", popt[0], "+/-", perr[0], " x^2 = ", rchi_sq
    
    # drop several points
    len_x = len(indat[0])
    idrop=int(len_x*0.1)
    xdata[ifile] = indat[0][idrop:int(len_x*0.88)] # also drop some tail points
    ydata = np.zeros(len(xdata[ifile]))
    for i in range(len(xdata[ifile])):
        if np.abs(indat[1][idrop+i]) > 1e-13:
            ydata[i] = np.log(np.abs(indat[1][idrop+i]))
        else :
            break
    ndata = i+1
    xdata[ifile] = xdata[ifile][range(ndata)]
    ydata = ydata[range(ndata)]
    ydata_sigma = len(xdata[ifile])*[1]
    #p0=2*[1.0] #
    #popt,perr,rchi_sq = chi_square( xdata, ydata, ydata_sigma, p0 )
    #print "drop 1point, fitting results :"
    #print "a = ", popt[0], "+/-", perr[0]
    #print "b = ", popt[1], "+/-", perr[1]
    #print " x^2 = ", rchi_sq
    
    p0=2*[1.0] #
    popt[ifile],perr,rchi_sq = chi_square( xdata[ifile], ydata, ydata_sigma, p0 )
    #print popt[0], "+",popt[1],"*log(x)"
    tag=(" %7.3f *x**( %7.3f )" % (np.exp(popt[ifile][0]), popt[ifile][1]) )
    print tag.replace(" ", "")
    #print " x^2 = ", rchi_sq
    
    # normalization
    ydata_fit = curvefunc(xdata[ifile], *popt[ifile])
    ydata_norm[ifile] = indat[1][idrop:ndata+idrop]*np.exp(-ydata_fit)
    
    # fft
    nk = 200
    sp[ifile] = np.fft.fft(ydata_norm[ifile], n=nk) / nk
    freq[ifile] = np.fft.fftfreq(nk)
    sp[ifile] = sp[ifile][range(nk/2)]
    freq[ifile] = freq[ifile][range(nk/2)]
    
    #orignial data
    x1[ifile]=indat[0]
    y1[ifile]=indat[1]

    # normalized data
    pos_y[ifile] = ydata_norm[ifile].copy()
    neg_y[ifile] = ydata_norm[ifile].copy()
    pos_y[ifile][pos_y[ifile] <= 0] = np.nan
    neg_y[ifile][neg_y[ifile] > 0] = np.nan
    
fig = plt.figure(figsize=(11, 8)) 
gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1]) 


mcolor=['b','g','r', 'c','m','y']
#mcolor=['b','g','r', 'g','b','r']
mmarker=['o','s','d', 'v', '^', '<']
flabel=["$,\ \\frac{1}{8}$ doping", "$,\ \\frac{1}{12}$ doping", "$,\ \\frac{1}{16}$ doping",
        "$,\ \\frac{1}{8}$ doping", "$,\ \\frac{1}{12}$ doping", "$,\ \\frac{1}{16}$ doping", ]

ax0 = plt.subplot(gs[0,0])
ax0.set_xlabel('$x-x_0$', fontsize=14)
ax0.set_ylabel('$|S(x-x_0)|$', fontsize=14)
for ifile in range(nfile/2):
    ax0.loglog(x1[ifile],y1[ifile]*10**(-2.0*ifile), mmarker[ifile], markersize=8, markerfacecolor=mcolor[ifile], markeredgecolor=mcolor[ifile])
    ax0.loglog(x1[ifile],-y1[ifile]*10**(-2.0*ifile), mmarker[ifile], markersize=8, markerfacecolor='none', markeredgecolor=mcolor[ifile])
    x0=np.arange(1.0,len(x1[ifile]),0.01)
    ax0.loglog(x0,np.exp(popt[ifile][0])*x0**popt[ifile][1]*10**(-2.0*ifile), '--', color=mcolor[ifile],
    label="$f(x-x_0)=%7.3f$"%np.exp(popt[ifile][0])+"$(x-x_0)$"+"$^{%7.3f}$"%popt[ifile][1]+flabel[ifile])

ax0.set_xticks(np.arange(5,40,5))
minor_ticks=np.arange(10)/20.0
ax0.set_xticks(minor_ticks, minor=True)
ax0.get_xaxis().set_major_formatter(ScalarFormatter())

ax0.set_xlim([idrop-0.5,int(len_x*0.88)-0.5])
#ax0.set_xlim([1,len_x])
ax0.set_ylim([0.1**12,1.0])
#ax[0].set_yscale('log')
ax0.grid(which='major')
#ax0.grid(which='minor')
ax0.text(-0.12, 1.02, 'a', transform=ax0.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax0.legend(fontsize=12)

ax3 = plt.subplot(gs[0,1])
ax3.set_xlabel('$x-x_0$', fontsize=14)
ax3.set_ylabel('$|D(x-x_0)|$', fontsize=14)
for ifile in range(nfile/2,nfile):
    ax3.loglog(x1[ifile],y1[ifile]*10**(-2.0*(ifile-3)), mmarker[ifile], markersize=8, markerfacecolor=mcolor[ifile], markeredgecolor=mcolor[ifile])
    ax3.loglog(x1[ifile],-y1[ifile]*10**(-2.0*(ifile-3)), mmarker[ifile], markersize=8, markerfacecolor='none', markeredgecolor=mcolor[ifile])
    x0=np.arange(1.0,len(x1[ifile]),0.01)
    ax3.loglog(x0,np.exp(popt[ifile][0])*x0**popt[ifile][1]*10**(-2.0*(ifile-3)), '--', color=mcolor[ifile],
    label="$f(x-x_0)=%7.3f$"%np.exp(popt[ifile][0])+"$(x-x_0)$"+"$^{%7.3f}$"%popt[ifile][1]+flabel[ifile])

ax3.set_xticks(np.arange(5,40,5))
minor_ticks=np.arange(10)/20.0
ax3.set_xticks(minor_ticks, minor=True)
ax3.get_xaxis().set_major_formatter(ScalarFormatter())

ax3.set_xlim([idrop-0.5,int(len_x*0.88)-0.5])
#ax3.set_xlim([1,len_x])
ax3.set_ylim([0.1**14,1.0])
#ax[0].set_yscale('log')
ax3.grid(which='major')
#ax3.grid(which='minor')
ax3.text(-0.12, 1.02, 'b', transform=ax3.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax3.legend(fontsize=12)

 
ax1 = plt.subplot(gs[1,0])
ax1.set_xlabel('$x-x_0$', fontsize=14)
ax1.set_ylabel('$S\ or\ D(x-x_0)/f(x-x_0)$',fontsize=14)
ax1.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
ax1.grid(which='major')
ax1.grid(which='minor')
ax1.text(-0.12, 1.02, 'c', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

for ifile in range(nfile):
    ax1.plot(xdata[ifile],ydata_norm[ifile]-5.0*ifile, '-', color=mcolor[ifile], linewidth=0.1)
    ax1.plot(xdata[ifile], pos_y[ifile]-5.0*ifile, mmarker[ifile], markersize=8, markerfacecolor=mcolor[ifile],   markeredgecolor=mcolor[ifile])
    ax1.plot(xdata[ifile], neg_y[ifile]-5.0*ifile, mmarker[ifile], markersize=8, markerfacecolor='none', markeredgecolor=mcolor[ifile])

ax2 = plt.subplot(gs[1,1])
ax2.set_xlabel('$qa_0/2\pi$',fontsize=14)
ax2.set_ylabel('|F(q)|',fontsize=14)
ax2.set_xticks(minor_ticks, minor=True)
ax2.grid(which='major')
ax2.grid(which='minor')
ax2.text(-0.12, 1.02, 'd', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

for ifile in range(nfile):
    ax2.plot(freq[ifile],abs(sp[ifile])-0.1*ifile,mcolor[ifile]) # plotting the spectrum

plt.tight_layout()
plt.savefig("fig3.pdf", dpi=300)
