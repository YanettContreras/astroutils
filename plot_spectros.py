import pyfits as pf
import numpy as n
import pylab as plb
import scipy.signal as sss
import gaussfit2d as g2d
import os
import time
import string as st
import gaussfitting as gft
from matplotlib.ticker import MaxNLocator
import matplotlib

def plotspectros(filein,xc,yc,x1,x2,y1,y2,dd=1, gf=True, data=False, datosx=0, aa=30):
    matplotlib.rc('xtick', labelsize=8)
    matplotlib.rc('ytick', labelsize=8)
    if data==True:
        fitsf=pf.getdata(filein)
        fitsh=pf.getheader(filein)
        nx1 = fitsh["NAXIS1"]
        delt1 = fitsh["CDELT1"]
        ra = delt1*(n.arange(nx1)+1-fitsh["CRPIX1"]) + fitsh["CRVAL1"]
        nx2 = fitsh["NAXIS2"]
        delt2 = fitsh["CDELT2"]
        dec = delt2*(n.arange(nx2)+1-fitsh["CRPIX2"]) + fitsh["CRVAL2"]
        nx3 = fitsh["NAXIS3"]
        delt3 = fitsh["CDELT3"]
        velo = (delt3*(n.arange(nx3)+1-fitsh["CRPIX3"]) + fitsh["CRVAL3"])/1000
    else:
        fitsf=filein
        velo=datosx
    mmx=n.max(fitsf[:,xc,yc])
    print 'Maximo en xc, yc :',mmx
    dx=x2-x1
    dy=y2-y1
    
    i=0
    for j in range(0,dx/dd)[::-1]:
        for k in range(0,dy/dd):
            i=i+1
            ax=plb.subplot(dx/dd,dy/dd,i)
            ax.xaxis.set_major_locator(MaxNLocator(3))
            ax.yaxis.set_major_locator(MaxNLocator(5))
            if i!=((dy/dd)*((dx/dd)-1)+1):
                ax.set_yticklabels([])
                ax.set_xticklabels([])
            #print (j*dd)+x1, (k*dd)+y1, i
            #plb.plot(velo, sss.medfilt(fitsf[:,(j*dd)+x1,(k*dd)+y1], kernel_size=3))
            plb.plot(velo, fitsf[:,(j*dd)+x1,(k*dd)+y1])
            plb.ylim(-.2,mmx*1.2)
            #print n.max(sss.medfilt(fitsf[:,(j*dd)+x1,(k*dd)+y1], kernel_size=3))
            if gf==True:
                if (j*dd)+x1==xc and (k*dd)+y1==yc:
                    #print n.max(fitsf[:,xc,yc])
                    #print n.max(fitsf[:,(j*dd)+x1,(k*dd)+y1])
                    xx,ww,mm=doplotgauss(velo, fitsf[:,(j*dd)+x1,(k*dd)+y1], ax, aa=aa)
                    plb.ylim(-.2,mmx*1.2)
    try:
        #plb.text(-100,-.5,'C: %5.2f W:%5.2f P:%5.2f'%(xx,ww,mm),horizontalalignment='center',verticalalignment='top')
        plb.show()
        
        return xx,ww,mm
    except UnboundLocalError:
        print "Grid no contiene al mamximo"
        return 0,0,0

  


def doplotgauss(datosx, datosy, axi, putval=False, aa=30):
    centrolam=sss.medfilt(datosy, kernel_size=5)
    center=n.argmax(centrolam)
    fit, xxx, ww, mm= gft.gaussfitm (datosx, centrolam , center,aa=aa)
    s1=gft.gaussfit (datosx, datosy , center)
    fit2 = lambda t : s1[2]*n.exp(-(t-s1[0])**2/(2*s1[1]**2))
    print "RESULTADOS GAUSSIANA ESPECTRO CENTRAL"
    print xxx,ww,mm, n.max(datosy)
    print s1
    #ax2.plot(velo[maxl-ll:maxl+ll], fit(velo[maxl-ll:maxl+ll]), 'k')
    if n.abs(s1[0]-xxx)>100 or n.abs(s1[1]-ww)>1:
        axi.plot(datosx, fit(datosx), 'r')
        if putval==True:
            axi.text(datosx[-1]-0.09*ww,0.7*mm,'C: %5.2f\nW:%5.2f\nP:%5.2f'%(xxx, ww, n.max(datosy)), fontsize='xx-small')
        return xxx, ww, n.max(datosy)
    else:
        axi.plot(datosx, fit2(datosx), 'r')
        if putval==True:
            axi.text(datosx[-1]-0.09*ww,0.7*mm,'C: %5.2f\nW:%5.2f\nP:%5.2f'%(s1[0], s1[1], n.max(datosy)), fontsize='xx-small')
        return s1[0], s1[1], n.max(datosy)


