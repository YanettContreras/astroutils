import numpy as n
from curva_draine import *
import scipy.stats as stt

def bbmod(Tem,nu):
    nu=nu*1e9
    c=2.997924e10 #cm/s
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    bb=(2*h*nu**3/(c**2))*(1/(n.exp(h*nu/(k*Tem))-1))
    #bb=2*k*Tem*nu**3/(c**2)
    return bb


def bbmod2(Tem,la):
    la=la*1e-4
    c=2.997924e10 #cm/s
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    bb=(2*h*(c**2)/(la**5))*(1/(n.exp(h*c/(la*k*Tem))-1))
    return bb


def Nh2(Speak,la,Tem,beam):
    """
    N(H2)_nu= S_nu(peak)/(beam*B(Tc)*sigma)
    From Garay, G. et at 2007
    """
    kla=curva_draine(la)
    Speak=Speak*1e-23  #Jy to erg/s/cm^2/Hz
    beam=n.pi*(beam/206265.)**2    #Pi r^2 en radians
    Rdg=1/125.  #Draine 2003
    mu=2.8     #Adopting 10% abundance of He with respect to H
    mh=1.67e-24   #gr
    sigma=Rdg*mu*mh*kla
    ff=lamdatofrec(la)
    return Speak/(beam*bbmod(Tem,ff)*sigma)

def Nh2870(Speak,la,Tem,beam):
    """
        N(H2)_nu= S_nu(peak)/(beam*B(Tc)*sigma)
        From Garay, G. et at 2007
        """
    kla=0.012
    Speak=Speak*1e-23  #Jy to erg/s/cm^2/Hz
    beam=n.pi*(beam/206265.)**2    #Pi r^2 en radians
    Rdg=1/125.  #Draine 2003
    mu=2.8     #Adopting 10% abundance of He with respect to H
    mh=1.67e-24   #gr
    sigma=Rdg*mu*mh*kla
    ff=lamdatofrec(la)
    return Speak/(beam*bbmod(Tem,ff)*sigma)


def Nh2new(Speak,la,Tem,beam):
    """
        N(H2)_nu= S_nu(peak)/(beam*B(Tc)*sigma)
        From Garay, G. et at 2007
        """
    kla=0.012
    Speak=Speak*1e-23  #Jy to erg/s/cm^2/Hz
    beam=np.pi/(4.*np.log(2)*(206265.**2.))*beam*3600*beam*3600
    #Pi r^2 en radians
    Rdg=1/125.  #Draine 2003
    mu=2.8     #Adopting 10% abundance of He with respect to H
    mh=1.67e-24   #gr
    sigma=Rdg*mu*mh*kla
    ff=lamdatofrec(la) #la in micron
    return Speak/(beam*bbmod(Tem,ff)*sigma)




def Fhot(Thot,Twar,Tcol,Rhot,Rw,Rc,la,Speak,beam):
    """
    Tem= T
    R: Radius of the source in arcsec
    """
    
    kla=curva_draine(la)
    ff=lamdatofrec(la)
    s_ang=n.pi*(Rhot/206265.)**2###n.pi*(R/D)**2
    s_angw=n.pi*(Rw/206265.)**2
    s_angc=n.pi*(Rc/206265.)**2
    Rdg=1/125.  #Draine 2003
    mu=2.8     #Adopting 10% abundance of He with respect to H
    mh=1.67e-24   #gr
    sigma=Rdg*mu*mh*kla
    Ncol=Nh2(Speak,1200.,Tcol,beam)
    Nwar=n.sqrt(s_angw/s_angc)*Ncol
    Nhot=n.sqrt(s_ang/s_angc)*Ncol
    fac1=(1-n.exp(-Nhot*sigma))
    fac2=n.exp(-(Ncol+Nwar)*sigma/2.)
    f_nu=s_ang*bbmod(Thot,ff)*fac1*fac2
    return f_nu

def Fwarm(Twar,Tcol,Rw,Rc,la,Speak,beam):
    kla=curva_draine(la)
    ff=lamdatofrec(la)
    R=Rw/206265.
    R2=Rc/206265.
    s_ang=n.pi*R**2###n.pi*(R/D)**2
    s_angc=n.pi*R2**2###n.pi*(R/D)**2
    Rdg=1/125.  #Draine 2003
    mu=2.8     #Adopting 10% abundance of He with respect to H
    mh=1.67e-24   #gr
    sigma=Rdg*mu*mh*kla
    Ncol=Nh2(Speak,1200.,Tcol,beam)
    Nwar=n.sqrt(s_ang/s_angc)*Ncol
    fac1=(1-n.exp(-Nwar*sigma))
    fac2=n.exp(-Ncol*sigma/2)
    f_nu=s_ang*bbmod(Twar,ff)*fac1*fac2
    return f_nu

def Fcol(Tcol,Rc,la,la0,beta1):
    """
    Tcol= Temperature cold component
    la  = Frecuencia en Hz
    beta1= exp of the gray body 
    Returns Fcol as F_nu  in cgs units (to convert to Jy divide it
    by 1e-23
    """
    ff=lamdatofrec(la)
    ff0=lamdatofrec(la0)
    R=Rc/206265.
    s_ang=n.pi*R**2###n.pi*(R/D)**2
    fac1=(1-n.exp(-(ff/ff0)**beta1))
    f_nu=s_ang*bbmod(Tcol,ff)*fac1
    return f_nu

def Ftot(dthot,dtwar,Tcol,drw,drc,Rhot,la0,beta1,la,Speak,beam):
    """
    Returns Ftot in cgs
    """
    Twar=dtwar+Tcol
    Thot=dthot+Twar
    Rw=Rhot+drw
    Rc=Rw+drc
    hot=Fhot(Thot,Twar,Tcol,Rhot,Rw,Rc,la,Speak,beam)
    warm=Fwarm(Twar,Tcol,Rw,Rc,la,Speak,beam)
    cold=Fcol(Tcol,Rc,la,la0,beta1)
    return hot+warm+cold

def Fcw(dtwar,Tcol,Rw,drc,la0,beta1,la,Speak,beam):
    """
    Returns Ftot in Jy
    """
    Twar=dtwar+Tcol
    Rc=Rw+drc
    warm=Fwarm(Twar,Tcol,Rw,Rc,la,Speak,beam)
    cold=Fcol(Tcol,Rc,la,la0,beta1)
    return warm+cold

def lamdatofrec(lamda):
    frec=(3.e8/(lamda*1.e-6))
    frec=frec/1.e9
    return frec

def frectolamda(nu):
    la=(3.e8/(nu*1.e9))
    return la

def modbb(Tem,la,la0,beta1,Rc):
    ff=lamdatofrec(la)
    ff0=lamdatofrec(la0)
    R=Rc/206265.
    s_ang=R**2###n.pi*(R/D)**2
    fac1=(ff/ff0)**beta1
    f_nu=s_ang*bbmod(Tem,ff)*fac1
    return f_nu


def modbb_l(Tem,la,la0,beta1,Rc):
    ff=lamdatofrec(la)
    ff0=lamdatofrec(la0)
    R=Rc/206265.
    s_ang=R**2###n.pi*(R/D)**2
    fac1=(ff/ff0)**beta1
    f_nu=s_ang*bbmod(Tem,ff)*fac1
    f_lambda=f_nu*2.9979e10/((la*1e-4)**2)
    return f_lambda


def sed2c(Tcol,Thot,beta1,beta2,Rc,la,la0):
    cold=modbb(Tcol,la,la0,beta1,Rc)
    hot=modbb(Thot,la,la0,beta2,Rc)
    return (cold+hot)/1.e-23

def fit_2c(args,Rc,la,la0,fluxes):
    ff=lamdatofrec(la)
    Tcol,Thot,beta1,beta2=args
    sed1=sed2c(Tcol,Thot,beta1,beta2,Rc,la,la0)
    chi2=stt.chisquare(fluxes*ff,sed1*ff)[0]
    return n.abs(chi2)
    
def fit_ftot(args,Rhot,la0,beta1,la,Speak,beam,fluxes):
    dthot,dtwar,Tcol,drw,drc=args
    tot1=Ftot(dthot,dtwar,Tcol,drw,drc,Rhot,la0,beta1,la,Speak,beam)
    ff=lamdatofrec(la)
    chi2=stt.chisquare(n.log(fluxes*ff),n.log(tot1*ff*1e23))[0]
    Twar=Tcol+dtwar
    Thot=Twar+dthot
    Rw=Rhot+drw
    Rc=Rw+drc
    if dtwar<0 or dthot<0:
        chi2=1000.
    if drw<0 or drc<0:
        chi2=1000.
    print chi2,Thot,Twar,Tcol,Rhot,Rw,Rc
    return n.abs(chi2)
