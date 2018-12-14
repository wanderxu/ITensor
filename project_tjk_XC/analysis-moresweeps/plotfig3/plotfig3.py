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
    #print indat

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
    
fig = plt.figure(figsize=(11, 10)) 
gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1]) 


mcolor=['r','g','b', 'y','m','c']
#mcolor=['b','g','r', 'g','b','r']
mmarker=['d','s','o', '<', '^', 'v']
#flabel=["$,\ \\frac{1}{16}$ doping", "$,\ \\frac{1}{12}$ doping", "$,\ \\frac{1}{8}$ doping",
#        "$,\ \\frac{1}{16}$ doping", "$,\ \\frac{1}{12}$ doping", "$,\ \\frac{1}{8}$ doping", ]
flabel=["$,\ $1/16 doping", "$,\ $1/12 doping", "$,\ $1/8 doping",
        "$,\ $1/16 doping", "$,\ $1/12 doping", "$,\ $1/8 doping", ]

ax0 = plt.subplot(gs[0,0])
ax0.set_xlabel('$x$', fontsize=14)
ax0.set_ylabel('$|S(x)|$', fontsize=14)
for ifile in range(nfile/2):
    ax0.loglog(x1[ifile],y1[ifile]*10**(-2.0*ifile), mmarker[ifile], markersize=8, markerfacecolor=mcolor[ifile], markeredgecolor=mcolor[ifile])
    ax0.loglog(x1[ifile],-y1[ifile]*10**(-2.0*ifile), mmarker[ifile], markersize=8, markerfacecolor='none', markeredgecolor=mcolor[ifile])
    x0=np.arange(1.0,len(x1[ifile]),0.01)
    ax0.loglog(x0,np.exp(popt[ifile][0])*x0**popt[ifile][1]*10**(-2.0*ifile), '--', color=mcolor[ifile],
    label="$f(x)=%7.3f$"%np.exp(popt[ifile][0])+"$x$"+"$^{%7.3f}$"%popt[ifile][1]+flabel[ifile])

ax0.set_xticks(np.arange(5,40,5))
minor_ticks=np.arange(10)/20.0
ax0.set_xticks(minor_ticks, minor=True)
ax0.get_xaxis().set_major_formatter(ScalarFormatter())

ax0.set_xlim([idrop-0.5,int(len_x*0.88)-0.5])
#ax0.set_xlim([1,len_x])
ax0.set_ylim([0.1**8,10**4])
#ax[0].set_yscale('log')
ax0.grid(which='major')
#ax0.grid(which='minor')
ax0.text(-0.12, 1.02, 'a', transform=ax0.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax0.legend(fontsize=14)

ax0.annotate('$\\times 10^{-4}$', fontsize=14, color='b',
            xy=(4.5, 0.1**5.8), xytext=(3.0, 0.1**7.0), 
            arrowprops=dict(arrowstyle="->", fc='b', ec='b') )
ax0.annotate('$\\times 10^{-2}$', fontsize=14, color='g',
            xy=(5.5, 0.1**3.6), xytext=(6.0, 0.1**5.0), 
            arrowprops=dict(arrowstyle="->", fc='g', ec='g') )

ax3 = plt.subplot(gs[0,1])
ax3.set_xlabel('$x$', fontsize=14)
ax3.set_ylabel('$|D(x)|$', fontsize=14)
for ifile in range(nfile/2,nfile):
    ax3.loglog(x1[ifile],y1[ifile]*10**(-2.0*(ifile-3)), mmarker[ifile], markersize=8, markerfacecolor=mcolor[ifile], markeredgecolor=mcolor[ifile])
    ax3.loglog(x1[ifile],-y1[ifile]*10**(-2.0*(ifile-3)), mmarker[ifile], markersize=8, markerfacecolor='none', markeredgecolor=mcolor[ifile])
    x0=np.arange(1.0,len(x1[ifile]),0.01)
    ax3.loglog(x0,np.exp(popt[ifile][0])*x0**popt[ifile][1]*10**(-2.0*(ifile-3)), '--', color=mcolor[ifile],
    label="$f(x)=%7.3f$"%np.exp(popt[ifile][0])+"$x$"+"$^{%7.3f}$"%popt[ifile][1]+flabel[ifile])

ax3.set_xticks(np.arange(5,40,5))
minor_ticks=np.arange(10)/20.0
ax3.set_xticks(minor_ticks, minor=True)
ax3.get_xaxis().set_major_formatter(ScalarFormatter())

ax3.set_xlim([idrop-0.5,int(len_x*0.88)-0.5])
#ax3.set_xlim([1,len_x])
ax3.set_ylim([0.1**10,10**4])
#ax[0].set_yscale('log')
ax3.grid(which='major')
#ax3.grid(which='minor')
ax3.text(-0.12, 1.02, 'b', transform=ax3.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax3.legend(fontsize=14)

ax3.annotate('$\\times 10^{-4}$', fontsize=14, color='c',
            xy=(4.5, 0.1**6.0), xytext=(3.0, 0.1**7.5), 
            arrowprops=dict(arrowstyle="->", fc='c', ec='c') )
ax3.annotate('$\\times 10^{-2}$', fontsize=14, color='m',
            xy=(5.3, 0.1**4.2), xytext=(6.0, 0.1**5.3), 
            arrowprops=dict(arrowstyle="->", fc='m', ec='m') )

 
ax1 = gridspec.GridSpecFromSubplotSpec(6, 1, subplot_spec=gs[1,0], hspace=0.07 )
#ax1 = plt.subplot(gs[1,0])
ax1arr = [ plt.subplot(ax1[i]) for i in range(nfile) ]
ax1arr[0].text(-0.12, 1.02, 'c', transform=ax1arr[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax1arr[nfile-1].set_xlabel('$x$', fontsize=14)
fig.text(.035, .25, '$S\ or\ D(x)/f(x)$', ha='center', va='center', rotation='vertical', fontsize=14)

# hide the spines in the middle
for ifile in range(nfile-1):
    ax1arr[ifile].spines['bottom'].set_visible(False)
    ax1arr[ifile].tick_params(bottom=False)
for ifile in range(1,nfile):
    ax1arr[ifile].spines['top'].set_visible(False)
ax1arr[0].xaxis.tick_top()
ax1arr[0].tick_params(labeltop=False)  # don't put tick labels at the top
ax1arr[nfile-1].xaxis.tick_bottom()

shift=8.0
for ifile in range(nfile):
    ax1arr[ifile].set_ylim([-ifile*shift-3.0,-ifile*shift+3.0])
    ax1arr[ifile].set_yticks([-ifile*shift-2, -ifile*shift, -ifile*shift+2])
    ax1arr[ifile].set_yticklabels( ["-2", "0", "2"] )
ifile=2
ax1arr[ifile].set_ylim([-ifile*shift-5.0,-ifile*shift+6.0])
ax1arr[ifile].set_yticks([-ifile*shift-3, -ifile*shift, -ifile*shift+3])
ax1arr[ifile].set_yticklabels( ["-3", "0", "3"] )
ifile=1
ax1arr[ifile].set_ylim([-ifile*shift-2.5,-ifile*shift+2.5])
ifile=3
ax1arr[ifile].set_ylim([-ifile*shift-3.5,-ifile*shift+3.5])

# From https://matplotlib.org/examples/pylab_examples/broken_axis.html
d = .005  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
for ifile in range(nfile-1):
    kwargs = dict(transform=ax1arr[ifile].transAxes, color='k', clip_on=False)
    ax1arr[ifile].plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax1arr[ifile].plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
for ifile in range(1, nfile):
    kwargs.update(transform=ax1arr[ifile].transAxes)  # switch to the bottom axes
    ax1arr[ifile].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax1arr[ifile].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

for ifile in range(nfile):
    ax1arr[ifile].axhline(y=-ifile*shift, color='k', linestyle='--', linewidth=0.5)
    ax1arr[ifile].grid(which='major')
    ax1arr[ifile].grid(which='minor')

for ifile in range(nfile):
    ax1arr[ifile].plot(xdata[ifile],ydata_norm[ifile]-shift*ifile, '-', color=mcolor[ifile], linewidth=0.1)
    ax1arr[ifile].plot(xdata[ifile], pos_y[ifile]-shift*ifile, mmarker[ifile], markersize=8, markerfacecolor=mcolor[ifile],   markeredgecolor=mcolor[ifile])
    ax1arr[ifile].plot(xdata[ifile], neg_y[ifile]-shift*ifile, mmarker[ifile], markersize=8, markerfacecolor='none', markeredgecolor=mcolor[ifile])

ax2 = plt.subplot(gs[1,1])
ax2.set_xlabel('$qa_0/2\pi$',fontsize=14)
ax2.set_ylabel('|F(q)|',fontsize=14)
ax2.set_xticks(minor_ticks, minor=True)
ax2.grid(which='major')
#ax2.grid(which='minor')
ax2.text(-0.12, 1.02, 'd', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

for ifile in range(nfile):
    ax2.plot(freq[ifile],abs(sp[ifile])-0.1*ifile,mcolor[ifile]) # plotting the spectrum

shifttag=[ "$-0.1$", "$-0.2$", "$-0.3$", "$-0.4$", "$-0.5$"]
for ifile in range(1,nfile):
    ax2.annotate(shifttag[ifile-1], fontsize=14, color=mcolor[ifile],
                xy=(0.07, -0.1*ifile+0.01), xytext=(0.12, -0.1*ifile+0.04), 
                arrowprops=dict(arrowstyle="->", fc=mcolor[ifile], ec=mcolor[ifile]) )

plt.tight_layout()
plt.savefig("fig3.pdf", dpi=300)
