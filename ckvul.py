import numpy as np
from numpy import sqrt, pi, exp, linspace
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.optimize import leastsq
import math

np_gam2 = 81
np_u = 146
lg_gam2_min = -6
lg_u_min = -16
step = 0.2
lg_gam2_max = lg_gam2_min + (np_gam2-1)*step
lg_u_max = lg_u_min + (np_u-1)*step
gff = np.genfromtxt("gauntff3.dat",delimiter="  ")

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
    for n in range(nr):
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

data = np.loadtxt('ckvul.dat')
x = data[:,0]
y = data[:,1]
ye = data[:,2]
#give diameter of the source
diam = math.sqrt(0.18*0.14)
o = 2*math.pi*(1-math.cos(diam*math.pi/(360*3600)))
nr = len(x)

a=60000.

def flux(x, b,eta):
    return o*(1.0-eta*eta)*a*x*x*3.072357e7*(-np.expm1(-b*0.05436*np.power(a,-1.5)*np.power(x,-2.0)*gauntff(x,a)))

popt, pcov = curve_fit(flux, x, y, bounds=([1.e7,0.],[3.e9,1.]), sigma=ye, method='trf', diff_step=1e-4)
popt
perr = np.sqrt(np.diag(pcov))
perr

print (popt)
print (perr)

print (y)
print (flux(x, popt[0], popt[1]))

chi = 0.0

chi = np.sqrt(np.power((y - flux(x, popt[0], popt[1])),2.)/ye/ye)
print (np.sum(chi)/(len(x) - 2))
