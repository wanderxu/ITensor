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
import matplotlib.patches as patches

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

nfile=2
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
    
fig = plt.figure(figsize=(8, 8)) 
gs = gridspec.GridSpec(2, 2, height_ratios=[1.2, 1]) 

ax0 = plt.subplot(gs[0,:])
ax0.set_xlabel('$x$', fontsize=14)
ax0.set_ylabel('$|P(x)|$', fontsize=14)

#mcolor=['b','orange']
#mmarker=['o','D']
mcolor=['r','olive']
mmarker=['d','D']
flabel=[", fitting $|P_{22}|$", ", fitting $|P_{12}|$"]
for ifile in range(nfile):
    ax0.loglog(x1[ifile],y1[ifile]*10**(-2*ifile), mmarker[ifile], markersize=8, markerfacecolor=mcolor[ifile], markeredgecolor=mcolor[ifile])
    ax0.loglog(x1[ifile],-y1[ifile]*10**(-2*ifile), mmarker[ifile], markersize=8, markerfacecolor='none', markeredgecolor=mcolor[ifile])
    x0=np.arange(1.0,len(x1[ifile]),0.01)
    ax0.loglog(x0,np.exp(popt[ifile][0])*x0**popt[ifile][1]*10**(-2*ifile), '--', color=mcolor[ifile],
    label="$f(x)=%7.3f$"%np.exp(popt[ifile][0])+"$x$"+"$^{%7.3f}$"%popt[ifile][1]+flabel[ifile])

ax0.set_xticks(np.arange(5,40,5))
minor_ticks=np.arange(10)/20.0
ax0.set_xticks(minor_ticks, minor=True)
ax0.get_xaxis().set_major_formatter(ScalarFormatter())

ax0.set_xlim([idrop-0.5,int(len_x*0.88)-0.5])
#ax0.set_xlim([1,len_x])
ax0.set_ylim([0.1**14,0.1**2])
#ax[0].set_yscale('log')
ax0.grid(which='major')
ax0.grid(which='minor')
ax0.text(-0.075, 1.05, 'a', transform=ax0.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax0.legend(fontsize=12)
#ax0.add_patch( patches.Rectangle( (0.1, 0.1), 0.5, 0.5, transform=ax0.transAxes, fill=True, facecolor='y', edgecolor='y') ) 
ax0.annotate('$\\times 10^{-2}$', fontsize=14, color='olive',
            xy=(14.0, 0.1**9.5), xytext=(11.5, 0.1**11.5), 
            arrowprops=dict(arrowstyle="->", fc='olive', ec='olive') )
 
ax1 = plt.subplot(gs[1,0])
ax1.set_xlabel('$x$', fontsize=14)
ax1.set_ylabel('$P(x)/f(x)$',fontsize=14)
ax1.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
ax1.grid(which='major')
ax1.grid(which='minor')
ax1.text(-0.15, 1.05, 'b', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

for ifile in range(nfile):
    ax1.plot(xdata[ifile],ydata_norm[ifile], '-', color=mcolor[ifile], linewidth=0.1)
    ax1.plot(xdata[ifile], pos_y[ifile], mmarker[ifile], markersize=8, markerfacecolor=mcolor[ifile],   markeredgecolor=mcolor[ifile])
    ax1.plot(xdata[ifile], neg_y[ifile], mmarker[ifile], markersize=8, markerfacecolor='none', markeredgecolor=mcolor[ifile])

ax2 = plt.subplot(gs[1,1])
ax2.set_xlabel('$qa_0/2\pi$',fontsize=14)
ax2.set_ylabel('|F(q)|',fontsize=14)
ax2.set_xticks(minor_ticks, minor=True)
ax2.grid(which='major')
ax2.grid(which='minor')
ax2.text(-0.15, 1.05, 'c', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax2.set_ylim([-0.01,0.27])

for ifile in range(nfile):
    ax2.plot(freq[ifile],abs(sp[ifile]),mcolor[ifile]) # plotting the spectrum

plt.tight_layout()
plt.savefig("fig1.pdf", dpi=300)