import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib 
import math
import matplotlib.axes as ax
from matplotlib.ticker import MaxNLocator
from numpy import sqrt, pi, exp, linspace
from scipy.optimize import curve_fit
from decimal import *

np_gam2 = 81
np_u = 146
lg_gam2_min = -6
lg_u_min = -16
step = 0.2
lg_gam2_max = lg_gam2_min + (np_gam2-1)*step
lg_u_max = lg_u_min + (np_u-1)*step
gff = np.genfromtxt("gauntff3.dat",delimiter="  ")

def lagrange2(xg, yg, lg_gam2):
    yval = 0.
    for i in range(0, 4):
        l = 1.0
        for j in range(0, 4):
            if i != j:
                l = l * (lg_gam2 - xg[j])/(xg[i]-xg[j])
        yval = yval + yg[i] * l
    return yval

def gauntff2(x, a):
    for i in range(0,5):
        gf=np.full(i,1.)
        xg=np.full(i,1.)
        g=np.full(i,1.)
    lg_gam2=np.log10(1.57887e5/a)
    lg_u = np.log10(x*4.79925e-2/a)
    ipg2 = (np.minimum(np.maximum((lg_gam2+6)*5,1)-1,77)).astype(int)
    ipu = (np.minimum(np.maximum(((lg_u+(16))*5),1)-1,144)).astype(int)
    for j in range(0,4):
        xg[j] = (ipg2+j)*0.2 - 6.0
    for j in range(0,4):
        for k in range(0,4):
            g[k] = gff[ipu+j][ipg2+k]
        gf[j] = lagrange2(xg, g, lg_gam2)
    for j in range(0,4):
        xg[j] = (ipu+j)*0.2 - 16.0
    gaunt = lagrange2(xg, gf, lg_u)
    return gaunt

def lagrange(xg, yg, lg_gam2, m):
    yval = 0.
    for i in range(0, 4):
        l = 1.0
        for j in range(0, 4):
            if i != j:
                l = l * (lg_gam2[m] - xg[m][j])/(xg[m][i]-xg[m][j])
        yval = yval + yg[m][i] * l
    return yval

def gauntff(x, a):
#    nr = range(len(x))
    for i in range(nr):
        gaunt=np.full(nr,1.)
        one=np.full(nr,1.)
        gf=np.full((nr,4),1.)
        xg=np.full((nr,4),1.)
        g=np.full((nr,4),1.)
        lg_gam2=np.full(nr,np.log10(1.57887e5/a))
    lg_u = np.log10(x*4.79925e-2/a)
    ipg2 = (np.minimum(np.maximum((lg_gam2+(6*one))*5,1)-one,77*one)).astype(int)
    ipu = (np.minimum(np.maximum(((lg_u+(16*one))*5),1)-one,144*one)).astype(int)
    for i in range(0,nr):
        for j in range(0,4):
            xg[i][j] = (ipg2[i]+j)*0.2 - 6.0
    for i in range(0,nr):
        for j in range(0,4):
            for k in range(0,4):
                g[i][k] = gff[ipu[i]+j][ipg2[i]+k]
            gf[i][j] = lagrange(xg, g, lg_gam2, i)
    for i in range(0,nr):
        for j in range(0,4):
            xg[i][j] = (ipu[i]+j)*0.2 - 16.0
    for l in range(0,nr):
        gaunt[l] = lagrange(xg, gf, lg_u, l)
    return gaunt

plt.figure(0,figsize=(4.,4.))
matplotlib.rcParams.update({'font.size': 10})


gs2 = gridspec.GridSpec(3, 3)
gs2.update(left=0.16, right=0.96, bottom=0.16, top=0.96, wspace=0.05, hspace=0.0)
ax3 = plt.subplot(gs2[:-1, :])
data = np.loadtxt("ckvul.dat")
col0 = data[: , 0]
col1 = data[: , 1]
col2 = data[: , 2]
diam = 1.44*np.sqrt(0.18*0.14)
o = 2*math.pi*(1-math.cos(diam*math.pi/(360*3600)))
a = 1.01331331e+04
b = 1.87820971e+07
x = np.linspace(1.,100.,1000)
nr = len(x)
p = b*0.054361255*np.power(a,-1.5)*np.power(x,-2.0)*gauntff(x,a)
y = -o*a*x*x*3.072357e7*np.expm1(-p)
nr = len(col1)
p2 = b*0.054361255*np.power(a,-1.5)*np.power(col0,-2.0)*gauntff(col0,a)
difference = col1 + o*a*col0*col0*3.072357e7*np.expm1(-p2)
plt.plot(x,y,color='red')
plt.errorbar(col0, col1, xerr=0, yerr=col2, linestyle=' ')
ax3.set_xlim([1.,100.])
ax3.set_xticks([])
ax3.set_ylim([0.2,2.])
ax3.set_ylabel('Flux (mJy)')
ax3.yaxis.set_label_coords(-0.15,0.5)
plt.yscale('log')
plt.xscale('log')
ax4 = plt.subplot(gs2[-1, :])
ax4.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))
ax4.set_ylabel('res. Flux (mJy)')
ax4.set_xlim([1.,100.])
plt.xscale('log')
ax4.set_xlabel('Frequency (GHz)')
ax4.yaxis.set_label_coords(-0.15,0.5)
plt.errorbar(col0, difference, xerr=0, yerr=col2, linestyle=' ')
plt.axhline(y=0.,color='red')

plt.show()
