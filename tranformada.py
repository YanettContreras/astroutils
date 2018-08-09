import numpy as n
import pyfits as pf
import numpy.fft as nfft
import pylab as pl
import scipy.fftpack as sfft 
import os

def ordena(fits2):
    x=n.shape(fits2)[0]/2
    y=n.shape(fits2)[1]/2
    bloq1=n.rot90(n.rot90(fits2[0:x,0:y]))
    bloq2=n.rot90(n.rot90(fits2[x:,0:y]))
    bloq3=n.rot90(n.rot90(fits2[x:,y:]))
    bloq4=n.rot90(n.rot90(fits2[0:x,y:]))
    ordena=n.zeros(n.shape(fits2))
    ordena[0:x,0:y]=bloq1
    ordena[x:,0:y]=bloq2
    ordena[x:,y:]=bloq3
    ordena[0:x,y:]=bloq4
    return ordena

#namefits='ATLASGAL-18.0-snr.fits'
#corte=20

def crea_imagen(namefits,corte,nameout):#(namefits,corte,nameout='None'):
    print nameout
    if nameout=='None':
        print 'ouch'
        nameout=namefits[0:]
    
    print nameout
    fits=pf.getdata(namefits)
    fits4=n.zeros(n.shape(fits))
    #fits1=nfft.fft2(fits[300:700,53:-200])
    fits=n.where(n.isnan(fits),0,fits)
    fits1=n.fft.fft2(fits)#[1:,1:])
    fits1_real=ordena(fits1.real)
    fits1_ima=ordena(fits1.imag)
    y,x = n.indices(n.shape(fits1_real))
    razon=n.float(n.shape(fits1_real)[0])/n.shape(fits1_real)[1]
    print razon
    r= n.sqrt((x-n.shape(fits1_real)[1]/2)**2+((y-n.shape(fits1_real)[0]/2)/razon)**2)
    
    fits1_rmask=fits1_real*n.where(r>corte,1,0)
    fits1_imask=fits1_ima*n.where(r>corte,1,0)
    fits1.real=ordena(fits1_rmask)
    fits1.imag=ordena(fits1_imask)
    fits3=n.fft.ifft2(fits1)
    fits4=n.sqrt(fits3.real**2+fits3.imag**2)
    
    h=pf.getheader(namefits)
    os.system('rm -f test'+nameout)
    os.system('rm -f real'+nameout)
    os.system('rm -f imag'+nameout)
    os.system('rm -f pow'+nameout)
    pf.writeto('test'+nameout,fits4,header=h)
    pf.writeto('real'+nameout,fits3.real,header=h)
    pf.writeto('imag'+nameout,fits3.imag,header=h)
    pf.writeto('pow'+nameout,n.sqrt(fits1.real**2+fits1.imag**2),header=h)
    #pl.imshow(fits3.real,vmin=0,vmax=10)
    return fits4


def crea_imagen_tampering(namefits,corte,nameout='None'):
    if nameout=='None':
        nameout=namefits[0:]
    fits=pf.getdata(namefits)
    fits4=n.zeros(n.shape(fits))
    #fits1=nfft.fft2(fits[300:700,53:-200])
    fits=n.where(n.isnan(fits),0,fits)
    fits1=n.fft.fft2(fits)#[1:,1:])
    fits1_real=ordena(fits1.real)
    fits1_ima=ordena(fits1.imag)
    y,x = n.indices(n.shape(fits1_real))
    razon=n.float(n.shape(fits1_real)[0])/n.shape(fits1_real)[1]
    print razon
    r= n.sqrt((x-n.shape(fits1_real)[1]/2)**2+((y-n.shape(fits1_real)[0]/2)/razon)**2)
    
    fits1_rmask=fits1_real*n.where(r<corte,1,0)
    fits1_imask=fits1_ima*n.where(r<corte,1,0)
    fits1.real=ordena(fits1_rmask)
    fits1.imag=ordena(fits1_imask)
    fits3=n.fft.ifft2(fits1)
    fits4=n.sqrt(fits3.real**2+fits3.imag**2)
    
    h=pf.getheader(namefits)
    os.system('rm -f test'+nameout)
    os.system('rm -f real'+nameout)
    os.system('rm -f imag'+nameout)
    os.system('rm -f pow'+nameout)
    pf.writeto('test'+nameout,fits4,header=h)
    pf.writeto('real'+nameout,fits3.real,header=h)
    pf.writeto('imag'+nameout,fits3.imag,header=h)
    pf.writeto('pow'+nameout,n.sqrt(fits1.real**2+fits1.imag**2),header=h)
    #pl.imshow(fits3.real,vmin=0,vmax=10)
    return n.std(fits3.real)
