import numpy as n
import pylab as pl
import pyfits as pf
from scipy import optimize
from numpy import *
from mpfit import mpfit
import time
import os

def gfit2d(filein,peak=0, x0=5, y0=5, x1=-5, y1=-5):
    #pl.clf()
    fitsf=pf.getdata(filein)
    
    if len(n.shape(fitsf))>2:
        fitsf=fitsf[0]
    
    
    fitsh=pf.getheader(filein)
    fitsf=n.where(n.isnan(fitsf),0,fitsf)
    max1=n.max(fitsf)
    
    if x0==None:
        x0=int(raw_input("Entre x0: "))
        y0=int(raw_input("Entre y0: "))
        x1=int(raw_input("Entre x1: "))
        y2=int(raw_input("Entre y1: "))
       
    data=fitsf[x0:x1,y0:y1]
    #pl.imshow(data,vmin=0, vmax=max1)
    params = fittwodgaussian(data, peak)
    #print 'resultado fit',params
    fit = twodgaussian(params[:-2],params[-2],params[-1])
    gfit=fit(*indices(data.shape))

    nx1 = fitsh["NAXIS1"]
    delt1 = fitsh["CDELT1"]
    ra = delt1*(n.arange(nx1)+1-fitsh["CRPIX1"]) + fitsh["CRVAL1"]
    nx2 = fitsh["NAXIS2"]
    delt2 = fitsh["CDELT2"]
    dec = delt2*(n.arange(nx2)+1-fitsh["CRPIX2"]) + fitsh["CRVAL2"]
    #nx3 = fitsh["NAXIS3"]
    #delt3 = fitsh["CDELT3"]
    #velo = (delt3*(n.arange(nx3)+1-fitsh["CRPIX3"]) + fitsh["CRVAL3"])/1000
    print "---------------------------------------------------------------"
    print "RESULTADOS"
    print filein
    print "---------------------------------------------------------------"
    print "Center   : ", ra[n.round(params[-2])+y0], dec[n.round(params[-1])+x0]
    print "Width x,y: ", params[0]*fitsh["CDELT2"]*3600, params[1]*fitsh["CDELT2"]*3600 
    print "Maximum  : ", params[2]
    print "Angle    : ", params[3]
    
    paramreal=n.array([ra[n.round(params[4])+y0], dec[n.round(params[5])+x0], params[0]*fitsh["CDELT2"]*3600, params[1]*fitsh["CDELT2"]*3600, params[2], params[3] ])
    params[4]=params[4]+y0
    params[5]=params[5]+x0
    return paramreal, params,gfit

def gauss2d(center_x,center_y,width_x,width_y,height):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y) 
    return lambda x,y: height*n.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)
   
def moments(data, peak=0):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    #print 'yyyyy',int(y),n.shape(data)
    if int(y)>=21:
        return 10,10,4,4,0.2
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    print n.where(data==height)
    if peak==1:
        a,b=n.where(data==height)
        x,y = a[0],b[0]
        print x,y
    return x, y, width_x, width_y,height
    
def fitgaussian2d(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gauss2d(*p)(*indices(data.shape)) -
                                    data)
    p, success = optimize.leastsq(errorfunction, params)
    print 'Parametros momento: ',params
    print 'Parametros ajuste: ',p
    p[-1]=params[-1]
    return p


def fittwodgaussian(data, peak=0):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data, peak)
    #print params
    params2=n.append(params[2:], 0)
    #print params2
    #errorfunction = lambda p: ravel(twodgaussian(p,cx,cy)(*indices(data.shape))-data)
    def errorfunction(p,cx,cy):
        return ravel(twodgaussian(p,cx,cy)(*indices(data.shape))-data)
    
    p, success = optimize.leastsq(errorfunction, params2, args=(params[0],params[1]))
    print 'Parametros momento: ',params
    print 'Parametros ajuste: ',p
    
    p[-2]=params2[-2]
    p=n.append(p,params[0])
    p=n.append(p,params[1])
    return p




def twodgaussian(param,center_x,center_y,shape=None):
    width_x,width_y,amplitude,rota=param
    amplitude = float(amplitude)
    center_x = float(center_x)
    center_y = float(center_y)
    rota = pi/180. * float(rota)
    rcen_x = center_x * n.cos(rota) - center_y * n.sin(rota)
    rcen_y = center_x * n.sin(rota) + center_y * n.cos(rota)
    def rotgauss(x,y):
        xp = x * n.cos(rota) - y * n.sin(rota)
        yp = x * n.sin(rota) + y * n.cos(rota)
        g = amplitude*n.exp(-(((rcen_x-xp)/width_x)**2+((rcen_y-yp)/width_y)**2)/2.)
        return g
    if shape is not None:
        return rotgauss(*n.indices(shape))
    else:
        return rotgauss

def savefitgauss(namefile, nmom0, par_real):
    t1=time.time()
    fecha=time.ctime(t1)
    datHandle = open(namefile, 'a')
    datHandle.write("%s\n"%nmom0)
    datHandle.write("%s\n"%fecha)
    datHandle.write("%f %f %f %f %f %f\n" %(par_real[0], par_real[1], par_real[2], par_real[3], par_real[4], par_real[5]))
    datHandle.close()
    return 0


def fitg2(nmom0, ssf, ssv=0, peak=0, vmax=2):
    pl.show()
    os.system('fits in='+nmom0+' out=tmpmom0.fits op=xyout')
    fmom0=pf.getdata('tmpmom0.fits')[0]
    pl.imshow(fmom0, vmin=0, vmax=1)
    par_real, param, gfit=gfit2d('tmpmom0.fits', x0=ssf, y0=ssf, x1=ssf, y1=ssf, peak=peak)
    #ssf=8
    pl.clf()
    pl.imshow(fmom0[ssf:-ssf,ssf:-ssf], vmin=0, vmax=vmax)
    pl.colorbar()
    pl.contour(gfit)
    if ssv==1:
        pl.savefig(nmom0[:-3]+'.ps')
        print "Archivo guardado: "+nmom0[:-3]+'.ps'
    return par_real
