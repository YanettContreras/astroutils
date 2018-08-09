import numpy as n
import scipy.special as ssp
import scipy.integrate as sitr
import scipy.interpolate as sit

def F_draine(a1,beta1,at):
    if beta1>=0:
        fd= 1+beta1*a1/at
    if beta1<0:
        fd= (1-beta1*a1/at)**(-1)
    return fd

def B_i(bc1=3e-5):
    """
    Usamos bc1=3e-5, para usar valor distinto
    usar la funcion como:
    B_i(bc1=new)
    """
    mc=1.994e-23
    rho1=2.24  #g cm^-3
    bci=n.array([0.75,0.25])*bc1
    ai=n.array([3.5,30])   #A
    sigma=0.4
    fac0=3/((2*n.pi)**(3/2))
    fac1=n.exp(-4.5*sigma**2)/rho1*ai**3*sigma
    fac2=bci*mc/(1+ssp.erf((3*sigma/n.sqrt(2))+n.log(ai/3.5)/(sigma*n.sqrt(2))))
    return fac0*fac1*fac2

def D_a(a1,bcn=3e-5):
    ai=n.array([3.5,30])
    sigma=0.4
    fac1=n.exp((-1/2)*(n.log(a1/ai)/sigma)**2)
    D_a=(B_i(bc1=bcn)/a1)*fac1
    return n.sum(D_a)

def nrg_c(a1,atg,acg,Cg,alpg,betag):
    a2=n.arange(3.5,a1,0.1)
    if a1<atg:
        fac0=1
    if a1>atg:
        fac0=n.exp(-((a2-atg)/acg)**3)
    D_ai=map(lambda x: D_a(x), a2)
    dnrg_c=D_ai+(Cg/a2)*(a2/atg)**alpg*F_draine(a2,betag,atg)*fac0
    NRG_c=sitr.trapz(a2,dnrg_c)
    return NRG_c

def nrg_s(a1,ats,acs,Cs,alps,betas):
    a2=n.arange(3.5,a1,0.1)
    if a1<ats:
        fac0=1
    if a1>ats:
        fac0=n.exp(-((a2-ats)/acs)**3)
    dnrg_s=(Cs/a2)*(a2/ats)**alps*F_draine(a2,betas,ats)*fac0
    NRG_s=sitr.trapz(a2,dnrg_s)
    return NRG_s

    

def curva_draine(la):
    R55=n.loadtxt('k_5.5A.txt')
    lamb1=R55[:,0]
    kabs=R55[:,4]
    fkabs=sit.UnivariateSpline(lamb1[::-1],kabs[::-1],s=0.0)
    return fkabs(la)
    
