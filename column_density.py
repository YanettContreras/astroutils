
import scipy.optimize as sop
import numpy as np
from uncertainties import unumpy

def col_dens_linmol(Tex, J, B, R,mu, tau, nu,nutb=True):
    """
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) From sanhueza 2012, and references therein
        Tex: Excitation temperature
        J=lower level rotational quantum number
        B= rotational constant
        R=relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        mu=permanent dipole moment
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        nutb=True mean use the case of optically thin assumption .
    """
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    Ej=h*B*J*(J+1)

    a1=(3*k)/(8*(pi**3)*B*(mu**2)*R)
    a2=(Tex+(h*B/(3*k)))/(J+1.)
    a3=np.exp(Ej/(k*Tex))/(1-np.exp(-h*nu/(k*Tex)))
    
    #print 'numeros  :',a1, (h*B/(3*k)), (J+1), Ej/k, -h*nu/(k)
    #a4=np.integral(tau,v1,v2)
    if nutb==True:
        a4=tau*(1/(jota(Tex,nu)-jota(Tbg,nu))) #where tau=integral (tb dv)
        print "N = %2.4g * (1/(J(Tex,nu)-%2.4g)[(Tex+%2.4g)/%2.4g] * [exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1,jota(Tbg,nu),(h*B/(3*k)),(J+1),Ej/k,-h*nu/(k))
    else:
        a4=tau #where tau=integral (tau dv)
        print "N = %2.4g * [(Tex+%2.4g)/%2.4g] * [exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1,(h*B/(3*k)),(J+1),Ej/k,-h*nu/(k))
    return a1*a2*a3*a4


def col_dens_linmol_err(Tex, J, B, R,mu, tau, nu,dtex,dtau,nutb=True):
    """
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) From sanhueza 2012, and references therein
        Tex: Excitation temperature
        J=lower level rotational quantum number
        B= rotational constant
        R=relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        mu=permanent dipole moment
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        nutb=True mean use the case of optically thin assumption .
        """
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    Ej=h*B*J*(J+1)
    
    Tex=unumpy.uarray(Tex,dtex)
    tau=unumpy.uarray(tau,dtau)
    
    a1=(3*k)/(8*(pi**3)*B*(mu**2)*R)
    a2=(Tex+(h*B/(3*k)))/(J+1.)
    a3=unumpy.exp(Ej/(k*Tex))/(1-unumpy.exp(-h*nu/(k*Tex)))
    
    #print 'numeros  :',a1, (h*B/(3*k)), (J+1), Ej/k, -h*nu/(k)
    #a4=np.integral(tau,v1,v2)
    if nutb==True:
        a4=tau*(1/(jota2(Tex,nu)-jota2(Tbg,nu))) #where tau=integral (tb dv)
        print "N = %2.4g * (1/(J(Tex,nu)-%2.4g)[(Tex+%2.4g)/%2.4g] * [exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1,jota(Tbg,nu),(h*B/(3*k)),(J+1),Ej/k,-h*nu/(k))
    else:
        a4=tau #where tau=integral (tau dv)
        print "N = %2.4g * [(Tex+%2.4g)/%2.4g] * [exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1,(h*B/(3*k)),(J+1),Ej/k,-h*nu/(k))
    return a1*a2*a3*a4





def col_dens_gralmol(Qrot, gu, aul, Ej, Tex, tau, nu, R=1, nutb=True):
    """
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) From sanhueza 2012, and references therein
        Qrot: Partition function
        gu: Statistical weight of the upper level
        aul: Eistein coefficient for spontaneous emission
        Tex: Excitation temperature
        R=relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        nutb=True mean use the case of optically thin assumption .
        """
    
    #EJ is EJ/k
    #Ej=0.00216
    R=1
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=np.exp(Ej/Tex)/(1-np.exp(-h*nu/(k*Tex)))
    print a1
    #print 'numeros  :',a1*a2, Ej, -h*nu/(k)
    #a4=np.integral(tau,v1,v2)
    if nutb==True:
        a4=tau*(1/(jota(Tex,nu)-jota(Tbg,nu))) #where tau=integral (tb dv)
        print a1#,a2,np.mean(a4),Ej,-h,nu,k
        print "N = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1*a2*np.mean(a4),Ej,-h*nu/(k))
        print "N = [%2.4g*(1/(J(Tex,nu)-%2.4g)*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1*a2,jota(Tbg,nu),Ej,-h*nu/(k))
        
    else:
        a4=tau #where tau=integral (tau dv)
        print "N = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1*a2,Ej,-h*nu/(k))
    
    return a1*a2*a3*a4


def col_dens_gralmol_err(Qrot, gu, aul, Ej, Tex, tau, nu,dtex,dtau, R=1, nutb=True):
    """
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) From sanhueza 2012, and references therein
        Qrot: Partition function
        gu: Statistical weight of the upper level
        aul: Eistein coefficient for spontaneous emission
        Tex: Excitation temperature
        R=relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        nutb=True mean use the case of optically thin assumption .
        """
    
    #EJ is EJ/k
    #Ej=0.00216
    R=1
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    Tex=unumpy.uarray(Tex,dtex)
    tau=unumpy.uarray(tau,dtau)
    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=unumpy.exp(Ej/Tex)/(1-unumpy.exp(-h*nu/(k*Tex)))
    
    if nutb==True:
        a4=tau*(1/(jota2(Tex,nu)-jota2(Tbg,nu))) #where tau=integral (tb dv)
    
    else:
        a4=tau #where tau=integral (tau dv)

    return a1*a2*a3*a4




def col_dens_gralmolLTE(Qrot, gu, aul, Eu, Tex, tau, nu, R=1, nutb=False):
    """
        Returns column density from Jones et al 2007, and references therein
        Qrot: Partition function
        gu: Statistical weight of the upper level
        aul: Eistein coefficient for spontaneous emission
        Tex: Excitation temperature
        R=relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        """
    
    #EJ is EJ/k
    #Ej=0.00216
    R=1
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=np.exp(Eu/Tex)
    print a1
    #print 'numeros  :',a1*a2, Ej, -h*nu/(k)
    #a4=np.integral(tau,v1,v2)
    if nutb==True:
        a4=tau*(1/(jota(Tex,nu)-jota(Tbg,nu))) #where tau=integral (tb dv)
        print a1#,a2,np.mean(a4),Ej,-h,nu,k
        print "N = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1*a2*np.mean(a4),Eu,-h*nu/(k))
        print "N = [%2.4g*(1/(J(Tex,nu)-%2.4g)*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1*a2,jota(Tbg,nu),Eu,-h*nu/(k))
    
    else:
        a4=tau #where tau=integral (tau dv)
        print "N = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1*a2,Eu,-h*nu/(k))
    
    return a1*a2*a3*a4


def col_dens_gralmolLTE_err(Qrot, gu, aul, Eu, Tex, tau, nu,dtex,dtau, R=1, nutb=False):
    """
        Returns column density from Jones et al 2007, and references therein
        Qrot: Partition function
        gu: Statistical weight of the upper level
        aul: Eistein coefficient for spontaneous emission
        Tex: Excitation temperature
        R=relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        """
    
    #EJ is EJ/k
    #Ej=0.00216
    R=1
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    Tex=unumpy.uarray(Tex,dtex)
    tau=unumpy.uarray(tau,dtau)
    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=unumpy.exp(Eu/Tex)
    if nutb==True:
        a4=tau*(1/(jota2(Tex,nu)-jota2(Tbg,nu))) #where tau=integral (tb dv)

    else:
        a4=tau #where tau=integral (tau dv)

    return a1*a2*a3*a4





def col_dens_hnco(Tex, tau):
    """
        Returns the column density (Garden et
        al 1991) From sanhueza 2012, and references therein
        Tex: Excitation temperature
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        Asumin R=1 which is diferent only for molecules with hyperfine structure
        """
    #EJ is EJ/k
    nu=87.925252e9
    Ej=6.32957
    R=1
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    Qrot=np.sqrt((pi*((Tex*k)**3))/((h**3)*918.417805e9*11.071010e9*10.910577e9))
    gu=9.0
    aul=8.78011e-6
    
    
    a1=(8*pi*(nu**3))/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=np.exp(Ej/Tex)/(1-np.exp(-(h*nu)/(k*Tex)))
    a4=1/(jota(Tex,nu)-jota(Tbg,nu))
    
    #print 'numeros  :',a1*a2, Ej, -h*nu/(k)
    #a4=np.integral(tau,v1,v2)
    a5=tau #where tau=integral (tau dv)
    #print Qrot
    #print "N[HNCO(404-303)] = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(a1*a2*a4,Ej,-h*nu/(k))
    
    print "N[HNCO(404-303)] = [%2.4g*(1/(J(Tex,nu)-%2.4g)*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2),jota(Tbg,nu),Ej,-h*nu/(k))
    return a1*a2*a3*a4*a5

def col_dens_hnco_err(Tex, tau, dtex,dtau):
    """
        Returns the column density (Garden et
        al 1991) From sanhueza 2012, and references therein
        Tex: Excitation temperature
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        Asumin R=1 which is diferent only for molecules with hyperfine structure
        """
    #EJ is EJ/k
    nu=87.925252e9
    Ej=6.32957
    R=1
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    Tex=unumpy.uarray(Tex,dtex)
    tau=unumpy.uarray(tau,dtau)

    
    Qrot=unumpy.sqrt((pi*((Tex*k)**3))/((h**3)*918.417805e9*11.071010e9*10.910577e9))
    gu=9.0
    aul=8.78011e-6
    
    
    a1=(8*pi*(nu**3))/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=unumpy.exp(Ej/Tex)/(1-unumpy.exp(-(h*nu)/(k*Tex)))
    a4=1/(jota2(Tex,nu)-jota2(Tbg,nu))
    
    a5=tau #where tau=integral (tau dv)
    
    return a1*a2*a3*a4*a5




def col_density_THA(gu,aul,nu,l,tex, ii):
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    
    B=nu/(2*(l+1))
    
    a1=4.4e-5
    a2=(nu**2)/(gu*aul*B)
    a3=np.exp((h*nu/(k*tex))*(1-(l/2)))
    print a2
    a4=tex*ii
    return a1*a2*a3*a4

def col_density_THA_err(gu,aul,nu,l,tex, ii,dtex,dtau):
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    
    Tex=unumpy.uarray(tex,dtex)
    tau=unumpy.uarray(ii,dtau)
    
    B=nu/(2*(l+1))
    
    a1=4.4e-5
    a2=(nu**2)/(gu*aul*B)
    a3=unumpy.exp((h*nu/(k*tex))*(1-(l/2)))
    print a2
    a4=tex*ii
    return a1*a2*a3*a4





def col_dens_cch(Tex, tau, optthin=False):
    """
        C2H(1-0)
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) From sanhueza 2012, and references therein
        Tex: Excitation temperature
        R, which is relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        
        Calculate tau from hyperfine
        
        attached file
        # # # #c2h
        u=2 #1
        l=1 #0
        Aul=0.153E-05
        nu=87.31690e9
        nughz=87.31690
        B0=43674.5177e6
        eu=4.193
        
        Tex=rows['dust_temperature']
        #    int_intensity=rows['c2h_ii']*2.0
        int_intensity=rows['c2h_tastar']*rows['c2h_fwhm']*1.064*2.0
        exp_part=math.exp( ((h*nu)/(k*Tex))*(1-(l/2)) )
        coldens=4.4e-5*const*exp_part*Tex*int_intensity

        
        """
    
    #EJ is EJ/k
    Ej=0.00216
    
    nu=87.316925e9
    R=5/12.
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    Qrot=((k*Tex)/(h*43.674518e9))+1/3.
    gu=5.0
    aul=1.52757e-6
    
    
    a1=(8*pi*(nu**3))/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=np.exp(Ej/Tex)/(1-np.exp(-h*nu/(k*Tex)))
    if optthin==True:
        a4=tau/(jota(Tex,nu)-jota(Tbg,nu)) #asuming optically thin, prob best to calculate tau form hyperfine
        print "N(C2H(1-0)) = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2/(jota(Tex,nu)-jota(Tbg,nu))),Ej,-h*nu/(k))
        print "N[C2H(1-0)] = [%2.4g*(1/(J(Tex,nu)-%2.4g)*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2),jota(Tbg,nu),Ej,-h*nu/(k))
    
    else:
        a4=tau #where tau=integral (tau dv)
        print "N(C2H(1-0)) = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2),Ej,-h*nu/(k))
    
    return a1*a2*a3*a4


def col_dens_cch_err(Tex, tau,dtex,dtau, optthin=False):
    """
        C2H(1-0)
        see col_dens_cch
        """
    #EJ is EJ/k
    Ej=0.00216
    nu=87.316925e9
    R=5/12.
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    Tex=unumpy.uarray(Tex,dtex)
    tau=unumpy.uarray(tau,dtau)
    
    Qrot=((k*Tex)/(h*43.674518e9))+1/3.
    gu=5.0
    aul=1.52757e-6
    
    a1=(8*pi*(nu**3))/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=unumpy.exp(Ej/Tex)/(1-unumpy.exp(-h*nu/(k*Tex)))
    if optthin==True:
        a4=tau/(jota2(Tex,nu)-jota2(Tbg,nu)) #asuming optically thin, prob best to calculate tau form hyperfine
    
    else:
        a4=tau #where tau=integral (tau dv)

    return a1*a2*a3*a4




def col_dens_cch_2(Tex, tau, optthin=False):
    """
        C2H(1-0)
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) From sanhueza 2012, and references therein
        Tex: Excitation temperature
        R, which is relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        
        Calculate tau from hyperfine
        
        # # # #c2h
        u=2 #1
        l=1 #0
        Aul=0.153E-05
        nu=87.31690e9
        nughz=87.31690
        B0=43674.5177e6
        eu=4.193
        
        Tex=rows['dust_temperature']
        #    int_intensity=rows['c2h_ii']*2.0
        int_intensity=rows['c2h_tastar']*rows['c2h_fwhm']*1.064*2.0
        exp_part=math.exp( ((h*nu)/(k*Tex))*(1-(l/2)) )
        coldens=4.4e-5*const*exp_part*Tex*int_intensity

        """
    
    #EJ is EJ/k
    Ej=0.00216
    nu=87.316925e9
    R=5/12.
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    Tbg=2.73
    
    Qrot=((k*Tex)/(h*43.674518e9))+1/3.
    gu=5.0
    aul=1.52757e-6
    
    
    #a1=(8*pi*nu**3)/((c**3)*R)
    #a2=Qrot/(gu*aul)
    #a3=np.exp(Ej/Tex)/(1-np.exp(-h*nu/(k*Tex)))
    
    B=43.674518e9
    J=0
    a1=3*k/(8*(pi**3)*B*(0.8e-18**2)*R)
    a2=(Tex+(h*B/(3.*k)))/(J+1)
    a3=np.exp(Ej/Tex)/(1-np.exp(-h*nu/(k*Tex)))
    #print 'met2',a1*a2*a3*a4*a5


    
    
    if optthin==True:
        a4=tau/(jota(Tex,nu)-jota(Tbg,nu)) #asuming optically thin, prob best to calculate tau form hyperfine
        print "N(C2H(1-0)) = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2/(jota(Tex,nu)-jota(Tbg,nu))),Ej,-h*nu/(k))
        print "N[C2H(1-0)] = [%2.4g*(1/(J(Tex,nu)-%2.4g)*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2),jota(Tbg,nu),Ej,-h*nu/(k))
    
    else:
        a4=tau #where tau=integral (tau dv)
        print "N(C2H(1-0)) = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2),Ej,-h*nu/(k))
    
    return a1*a2*a3*a4




def col_dens_hc3n(Tex, tau):
    """
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) From sanhueza 2012, and references therein
        Tex: Excitation temperature
        R, which is relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        """
    
    #EJ is EJ/k
    Ej=19.6484
    R=1
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    nu=91.199796e9
    Tbg=2.73
    Qrot=((k*Tex)/(h*4.5490586e9))+1/3.
    gu=21.0
    aul=58.13e-6
    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=np.exp(Ej/Tex)/(1-np.exp(-h*nu/(k*Tex)))
    a4=1/(jota(Tex,nu)-jota(Tbg,nu))
    
    #print 'numeros  :',a1*a2, Ej, -h*nu/(k)
    #a4=np.integral(tau,v1,v2)
    a5=tau #where tau=integral (tau dv)
    print "N(HC3N(10-9)) = [%2.4g*(1/(J(Tex,nu)-%2.4g)*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2),jota(Tbg,nu),Ej,-h*nu/(k))
    return a1*a2*a3*a4*a5

def col_dens_hc3n_err(Tex, tau,dtex,dtau):
    """
       """
    
    #EJ is EJ/k
    Ej=19.6484
    R=1
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    nu=91.199796e9
    Tbg=2.73
    
    Tex=unumpy.uarray(Tex,dtex)
    tau=unumpy.uarray(tau,dtau)

    Qrot=((k*Tex)/(h*4.5490586e9))+1/3.
    gu=21.0
    aul=58.13e-6
    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=unumpy.exp(Ej/Tex)/(1-unumpy.exp(-h*nu/(k*Tex)))
    a4=1/(jota2(Tex,nu)-jota2(Tbg,nu))
    
 
    a5=tau #where tau=integral (tau dv)
 
    return a1*a2*a3*a4*a5






def col_dens_hhcs(Tex, tau):
    """
        H2CS (313-212)
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) From sanhueza 2012, and references therein
        Tex: Excitation temperature
        R, which is relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        
        From CDS
        101477.8095  0.0014 -3.6338 3   12.5385 21  46509 303 3 1 3       2 1 2        H2CS
        mu=1.6491
        A = 1.16395e-20 nu^3 Sg mu_g^2/gup
        , where nu in MHz and mu in debyes (D)
        Partition Function Parameters
        A / MHz	291613.3
        B / MHz	17698.995
        C / MHz	16652.499
        
        From Paper Dae-Kim et. al (http://iopscience.iop.org/0067-0049/131/2/483/pdf/51376.web.pdf)
        Sij=2.67
        
        """
    
    #EJ is EJ/k
    
    R=1
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    
    Ej=12.5385
    nu=101.4778095e9
    Tbg=2.73
    
    #Qrot=np.sqrt((pi*((Tex*k)**3))/((h**3)*918.417805e9*11.071010e9*10.910577e9))
    #print 'qrot1',Qrot
    
    Qrot=np.sqrt((pi*((Tex*k)**3))/((h**3)*291.613e9*17.698e9*16.652e9))
    gu=21.0
    
    aul=1.16395e-20*((nu/1.e6)**3)*2.67*(((1.64)**2)/gu)
    
    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=np.exp(Ej/Tex)/(1-np.exp(-h*nu/(k*Tex)))
    a4=1/(jota(Tex,nu)-jota(Tbg,nu))
    
    #print 'numeros  :',a1,a2, Ej, -h*nu/(k)
    #print Qrot
    #a4=np.integral(tau,v1,v2)
    a5=tau #where tau=integral (tau dv)
    print a4
    print "N(H2CS(313-212)) = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2*a4),Ej,-h*nu/(k))
    print "N[H2CS(312-212)] = [%2.4g*(1/(J(Tex,nu)-%2.4g)*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2),jota(Tbg,nu),Ej,-h*nu/(k))
    return a1*a2*a3*a4*a5

def col_dens_htrcn(Tex, tau):
    """
        H13CN
        
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) From Pitann et all, and references therein
        Tex: Excitation temperature
        R, which is relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        
        From CDS
        
        86339.9215  0.0001 -2.5474 2    0.0000  9  28501 101 1           0            HC-13-N, v=0
        mu=2.9852
        
        A = 1.16395e-20 nu^3 Sg mu_g^2/gup
        , where nu in MHz and mu in debyes (D)
        
        Partition Function Parameters
       
        B / MHz	43170.127
       
        FRom http://iopscience.iop.org/0004-637X/575/1/250/fulltext/
        Sg=1.67
        
        """
    
    #EJ is EJ/k
    R=0.6
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s

    Ej=0.0
    nu=86.33992000000e9
    Tbg=2.73
    gu=9.0
    
    #Qrot=np.sqrt((pi*((Tex*k)**3))/((h**3)*918.417805e9*11.071010e9*10.910577e9))
    #print 'qrot1',Qrot
    
    Qrot=((k*Tex)/(h*43170.127e6))+(1/3.)  #if it were a lineal molec
    #Qrot=109.62  #For 75K
    
    
    aul=1.16395e-20*((nu/1.e6)**3)*1.67*(((2.9852)**2)/gu)
    #aul=2.2256e-5
    moldata=np.loadtxt('Molecular-data-info/h13cnp.txt', usecols=(0,4,5))
    Qrot=np.sum(moldata[:7,1]*np.exp(-moldata[:7,2]/Tex)) #up to 8th level
    print 'r',Qrot
    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=np.exp(Ej/Tex)/(1-np.exp(-h*nu/(k*Tex)))
    a4=1/(jota(Tex,nu)-jota(Tbg,nu))
    
    a5=tau #where tau=integral (tau dv)
    #print 'met1',a1*a2*a3*a4*a5
    
    #B=43170.127e6
    #J=0
    #a1=3*k/(8*(pi**3)*B*(2.9852e-18**2)*R)
    #a2=(Tex+(h*B/(3.*k)))/(J+1)
    #print 'met2',a1*a2*a3*a4*a5
    
    print "N(H13CN) = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2*a4),Ej,-h*nu/(k))
    print "N[H13CN] = [%2.4g*(1/(J(Tex,nu)-%2.4g)*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2),jota(Tbg,nu),Ej,-h*nu/(k))
    return a1*a2*a3*a4*a5

def col_dens_hcn(Tex, tau, thin=False):
    
    
    """


    HCN(1-0)
        
        ASUMING OPTICALLY THIN!!!!!!!
        
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991) 
        Tex: Excitation temperature
        R, which is relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        
        From CDS
        
        88631.6022  0.0001 -2.5140 2    0.0000  9  27501 101 1           0            HCN, v=0
        mu=2.9852
        
        A = 1.16395e-20 nu^3 Sg mu_g^2/gup
        , where nu in MHz and mu in debyes (D)
        
        Partition Function Parameters
        
        B / MHz	44315.976
        
        Tex=rows['dust_temperature']
        int_intensity=rows['hcn_ii']*2.0
        exp_part=math.exp( ((h*nu)/(k*Tex))*(1-(l/2)) )
        coldens=4.4e-5*const*exp_part*Tex*int_intensity
        
        # # # #hcn
        u=1
        l=0
        Aul=2.4075e-05
        nu=88.6316023e9
        nughz=88.6316023
        B0=44315.976e6
        eu=4.25

        
        """
    
    #EJ is EJ/k
    R=1#0.6
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    B=44315.976e6
    
    Ej=0.0
    nu=88.63160220e9
    Tbg=2.73
    gu=9.0
    sg=1.67 #
    

    Qrot=((k*Tex)/(h*B))+(1/3.)  #if it were a lineal molec
    
    aul=1.16395e-20*((nu/1.e6)**3)*sg*(((2.9852)**2)/gu)
    #aul=2.4075e-05
    print aul
    moldata=np.loadtxt('/Users/yanett/inv/MALT90/Models_caselli/Molecular-data-info/hcn.txt', usecols=(0,4,5))
    if type(Tex)==float:
        Qrot2=np.sum(moldata[:7,1]*np.exp(-moldata[:7,2]/Tex)) #up to 8th level
    else:
        j=0
        Qrot2=np.ones(len(Tex))
        for i in Tex:
            
            
            Qrot2[j]=np.sum(moldata[:7,1]*np.exp(-moldata[:7,2]/i))
            j=j+1

            
    print 'HCN : r',Qrot,Qrot2
    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=np.exp(Ej/Tex)/(1-np.exp(-h*nu/(k*Tex)))
    print 1/(jota(Tex,nu)-jota(Tbg,nu))
    if thin==False:
        a4=1#/(jota(Tex,nu)-jota(Tbg,nu))
    else:
        a4=1/(jota(Tex,nu)-jota(Tbg,nu))
    a5=tau #where tau=integral (tau dv)
    print 'met1',a1*a2*a3*a4*a5
    
    B=43170.127e6
    J=0
    a1=3*k/(8*(pi**3)*B*(2.9852e-18**2)*R)
    a2=(Tex+(h*B/(3.*k)))/(J+1)
    print 'met2',a1*a2*a3*a4*a5
    
    print "N(HCN) = [%2.4g*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2*a4),Ej,-h*nu/(k))
    print "N[HCN] = [%2.4g*(1/(J(Tex,nu)-%2.4g)*exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] * Int(Tau,dv)"%(np.mean(a1*a2),jota(Tbg,nu),Ej,-h*nu/(k))
    return a1*a2*a3*a4*a5


def col_dens_hcn_err(Tex, tau,dtex,dtau, thin=False):
    
    
    """
        
        
        HCN(1-0)
        
        ASUMING OPTICALLY THIN!!!!!!!
        
        Returns the column density of a linear , rigid roto molecule (Garden et
        al 1991)
        Tex: Excitation temperature
        R, which is relative intensity of the brightest hyperfine transition, it is 1 if the moel dont have hyperfine. N2H+=5/9. ; C2H=5/12.
        tau=integral of the opacity (\int tau dv)
        nu=frequeancy in Hz
        
        From CDS
        
        88631.6022  0.0001 -2.5140 2    0.0000  9  27501 101 1           0            HCN, v=0
        mu=2.9852
        
        A = 1.16395e-20 nu^3 Sg mu_g^2/gup
        , where nu in MHz and mu in debyes (D)
        
        Partition Function Parameters
        
        B / MHz	44315.976
        
        Tex=rows['dust_temperature']
        int_intensity=rows['hcn_ii']*2.0
        exp_part=math.exp( ((h*nu)/(k*Tex))*(1-(l/2)) )
        coldens=4.4e-5*const*exp_part*Tex*int_intensity
        
        # # # #hcn
        u=1
        l=0
        Aul=2.4075e-05
        nu=88.6316023e9
        nughz=88.6316023
        B0=44315.976e6
        eu=4.25
        
        
        """
    
    #EJ is EJ/k
    R=1#0.6
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    B=44315.976e6
    
    Ej=0.0
    nu=88.63160220e9
    Tbg=2.73
    gu=9.0
    sg=1.67 #
    
    Tex=unumpy.uarray(Tex,dtex)
    tau=unumpy.uarray(tau,dtau)
    
    Qrot=((k*Tex)/(h*B))+(1/3.)  #if it were a lineal molec
    
    aul=1.16395e-20*((nu/1.e6)**3)*sg*(((2.9852)**2)/gu)
    #aul=2.4075e-05
    print aul

    
    a1=(8*pi*nu**3)/((c**3)*R)
    a2=Qrot/(gu*aul)
    a3=unumpy.exp(Ej/Tex)/(1-unumpy.exp(-h*nu/(k*Tex)))
    print 1/(jota2(Tex,nu)-jota2(Tbg,nu))
    if thin==False:
        a4=1#/(jota(Tex,nu)-jota(Tbg,nu))
    else:
        a4=1/(jota2(Tex,nu)-jota2(Tbg,nu))
    a5=tau #where tau=integral (tau dv)
    print 'met1',a1*a2*a3*a4*a5
    
    B=43170.127e6
    J=0
    a1=3*k/(8*(pi**3)*B*(2.9852e-18**2)*R)
    a2=(Tex+(h*B/(3.*k)))/(J+1)
    print 'met2',a1*a2*a3*a4*a5
    
    
    return a1*a2*a3*a4*a5





def jota(t,nu):
    """
        Function J that goes into the relation between tau and Tmb.
        t=temperature
        nu=frequency in Hz
    """
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    a1=h*nu/(k)
    a2=1/(np.exp(h*nu/(k*t))-1)
    return a1*a2

def jota2(t,nu):
    """
        Function J that goes into the relation between tau and Tmb.
        t=temperature
        nu=frequency in Hz
        """
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    a1=h*nu/(k)
    a2=1/(unumpy.exp(h*nu/(k*t))-1)
    return a1*a2


def col_dens_sio(Tmb, B, mu,Tex,J, nu):
    """
        Returns the column density of SiO. From sanhueza 2012, and references therein
        Tmb: Main beam temperature integrated intensity (\int Tmb dv)
        J=lower level rotational quantum number
        B= rotational constant
        mu=permanent dipole moment
        Tex= excitation temperature
        nu=frequeancy in Hz
        """
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    
    R=1
    Tbg=2.73
    Ej=h*B*J*(J+1)
    
    a1=3*k/(8*(pi**3)*B*(mu**2)*R)
    a2=(Tex+(h*B/(3*k)))/(J+1)
    a3=np.exp(Ej/(k*Tex))/(1-np.exp(-h*nu/(k*Tex)))
    #print Ej,Ej/k*Tex,  (1-np.exp(-h*nu/(k*Tex))), a1, a2, a3
    a4=1/(jota(Tex,nu)-jota(Tbg,nu))
    a5=Tmb#.copy()
    
    
    #print 'numeros  :', a1, (h*B/(3*k)), (J+1), Ej/k, -h*nu/(k), jota(Tbg,nu)
    print "N(SiO) = %2.4g * [(Tex+%2.4g)/%2.4g] * [exp(%2.4g/Tex)/(1-exp(%2.4g/Tex)] x 1/[J(Tex)-%2.4g] * Int(Tau,dv)"%(a1,(h*B/(3*k)),(J+1),Ej/k,-h*nu/(k), jota(Tbg,nu))
    print "J(Tex)=%2.4g * 1/[exp(%2.4g/Tex)-1]"%(h*nu/(k),h*nu/(k))
    return a1*a2*a3*a4*a5

def col_dens_sio_err(Tmb, B, mu,Tex,J, nu, dtex,dtmb):
 
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    
    R=1
    Tbg=2.73
    Ej=h*B*J*(J+1)
    
    Tex=unumpy.uarray(Tex,dtex)
    Tmb=unumpy.uarray(Tmb,dtmb)

    
    a1=3*k/(8*(pi**3)*B*(mu**2)*R)
    a2=(Tex+(h*B/(3*k)))/(J+1)
    a3=unumpy.exp(Ej/(k*Tex))/(1-unumpy.exp(-h*nu/(k*Tex)))
    a4=1/(jota2(Tex,nu)-jota2(Tbg,nu))
    a5=Tmb#.copy()
    
    return a1*a2*a3*a4*a5



def col_dens_colow(II, Tkin):
    a1=2.4e14
    a2=Tkin+0.92
    a3=1./(1.-np.exp(-5.53/Tkin))
    a4=II
    return a1*a2*a3*a4




def r(ratio,tex,b1,b2,j1,j2,nu1,nu2):
    """
        return ratio between opacities of two isotopomers.
        ratio= ratio betwenn the isotopomenrs abundances
        Tex= excitation temperature
        b1=b2=Rotational constant of both molecules
        j1=j2=lower level rotational quantum number for both molecules
        nu1=nu2=frequencies of both molecules
        """
    c=2.98e10 #cm s-1
    pi=3.141516
    k=1.38*1e-16 #erg K-1
    h=6.626*1e-27 #erg s
    e2=h*b2*j2*(j2+1)
    e1=h*b1*j1*(j1+1)
    a1=((k*tex/(h*b2))+1/3.)/((k*tex/(h*b1))+1/3.)
    a2=np.exp(e2/(k*tex))/np.exp(e1/(k*tex))
    a3=(1-np.exp(-h*nu1/(k*tex)))/(1-np.exp(-h*nu2/(k*tex)))
    return ratio*a1*a2*a3

def tau_isotop(tau,r,tmb1,tmb2):
    """
        return equation to solve opacity between two isotopomers.
        r is the ratio between the opacities
        tmb1, tmb2=are the main beam temperature of both isotopomers
    """
    return (1-np.exp(-tau/r))/(1-np.exp(-tau))-(tmb2/tmb1)

def solve_tau(tmb1,tmb2,ratio,tex,b1,b2,j1,j2,nu1,nu2):
    """
        Returns the opacity from the data of two isotopomers (e.g. HCO+ and H13CO+)
        tmb1, tmb2: Main beam temperature integrated intensity (\int Tmb dv)
        j1.j2=lower level rotational quantum number
        tex= excitation temperature
        ratio=ratio between the abundaneces
        nu1,nu2=frequeancy in Hz
        """
    r1=r(ratio,tex,b1,b2,j1,j2,nu1,nu2)
    tau=.1
    tau2=sop.root(tau_isotop, tau, args=(r1,tmb1,tmb2))
    #print tau2.x
    tau1=sop.fsolve(tau_isotop, tau, args=(r1,tmb1,tmb2), maxfev=1000)
    
    if np.isnan(tmb1) or np.isnan(tmb2):
        return np.nan, np.nan, np.nan
    else:
        return tau1[0], tau2.x[0], tau1[0]/r1



def makehyperfine(p,mol,vx):
    p1,v0,dv,p4=p
    
    if mol=='hcn':
        D = [-7.07,0,4.84]
        r1,r2,r3=np.array([3,5,1])/9.   # Loughnane et al 2013
    if mol=='n2hp':
        D = [-8.20,0,5.74]
        #r1,r2,r3=[1.667,2.333,1]  #Keto & Rybicki 2010
        r1,r2,r3=np.array([5.,7.,3.])/15.  #CLASS Manual, Forveille et al.
    if mol=='c2h':
        D = [-40.17,0]
        r1,r2=np.array([0.4,0.2])/.6  ## Padovani et al 2009


    Z = (vx-v0)/dv   
    EZ = np.exp(-4*np.log(2)*Z**2)
    tmv=p4*r2*EZ
    tant2=(p1/p4)*(1.-np.exp(-p4*r2))

    W = (vx-v0-D[0])/dv
    EW = np.exp(-4*np.log(2)*W**2)
    t1v=p4*r1*EW
    tant1=(p1/p4)*(1.-np.exp(-p4*r1))

    if mol!='c2h':
        U = (vx-v0-D[2])/dv
        EU = np.exp(-4*np.log(2)*U**2)
        t3v=p4*r3*EU
        tant3=(p1/p4)*(1.-np.exp(-p4*r3))
    else:
        t3v=vx*0
        tant3=0.

    tau=t1v+tmv+t3v

    tant=(p1/p4)*(1.-np.exp(-tau))

    

    return tant, [tant1,tant2,tant3]


def makehyperfine_wf(p,dv,v0,mol,vx):
    #Makes hyperfine with fixed widths and positions
    p1,p4=p
    
    if mol=='hcn':
        D = [-7.07,0,4.84]
        r1,r2,r3=np.array([3,5,1])/9.   # Loughnane et al 2013
    if mol=='n2hp':
        D = [-8.20,0,5.74]
        #r1,r2,r3=[1.667,2.333,1]  #Keto & Rybicki 2010
        r1,r2,r3=np.array([5.,7.,3.])/15.
    if mol=='c2h':
        D = [-40.17,0]
        r1,r2=np.array([0.4,0.2])/.6   ## Padovani et al 2009
    
    
    Z = (vx-v0)/dv
    EZ = np.exp(-4*np.log(2)*Z**2)
    tmv=p4*r2*EZ
    tant2=(p1/p4)*(1.-np.exp(-p4*r2))
    
    
    W = (vx-v0-D[0])/dv
    EW = np.exp(-4*np.log(2)*W**2)
    t1v=p4*r1*EW
    tant1=(p1/p4)*(1.-np.exp(-p4*r1))

    if mol!='c2h':
        U = (vx-v0-D[2])/dv
        EU = np.exp(-4*np.log(2)*U**2)
        t3v=p4*r3*EU
        tant3=(p1/p4)*(1.-np.exp(-p4*r3))
    else:
        t3v=vx*0
        tant3=0.

    tau=t1v+tmv+t3v
    
    tant=(p1/p4)*(1.-np.exp(-tau))
    
    
    return tant, [tant1,tant2,tant3]



def hyper_chi2_amplitudes(p,dv,v0,mol,vx,amplitudes):
    #Return the chi2 between the model of the hyperfine and the amplitudes provided
    tant, ta=makehyperfine_wf(p,dv,v0,mol,vx)
    chi=np.sum(np.abs(np.array(ta)-np.array(amplitudes)))
    if (p[0]/amplitudes[1])<0.1 or (p[0]/amplitudes[1])>30 or p[1]<0:
        chi=100.
    return chi


def taudv_n2hp(t1,t2,t3,dv,v0,vx):
    #p=p1,v0,dv,p4   where p4 is the opacity
    p=[1.,0.1]   #Initial guess
    psolved=sop.fmin(hyper_chi2_amplitudes,p,args=(dv,v0,'n2hp',vx,[t1,t2,t3]))
    return psolved, hyper_chi2_amplitudes(psolved,dv,v0,'n2hp',vx,[t1,t2,t3])

def taudv_hcn(t1,t2,t3,dv,v0,vx):
    #p=p1,v0,dv,p4   where p4 is the opacity
    p=[1.,0.1]  #Initial guess
    psolved=sop.fmin(hyper_chi2_amplitudes,p,args=(dv,v0,'hcn',vx,[t1,t2,t3]))
    return psolved, hyper_chi2_amplitudes(psolved,dv,v0,'hcn',vx,[t1,t2,t3])


def taudv_c2h(t1,t2,t3,dv,v0,vx):
    #p=p1,v0,dv,p4   where p4 is the opacity
    p=[1.,0.1]  #Initial guess
    psolved=sop.fmin(hyper_chi2_amplitudes,p,args=(dv,v0,'c2h',vx,[t1,t2,t3]))
    return psolved, hyper_chi2_amplitudes(psolved,dv,v0,'c2h',vx,[t1,t2,t3])



def fillingfactor(tau,tex,tmb,nu):
    #print tau, tex, tmb, nu
    try:
        expfactor=np.exp(-tau)
    except:
        expfactor=np.array([np.exp(-x) for x in tau])
    ff=tmb/((jota(tex,nu)-jota(2.73,nu))*(1-expfactor))
    ff=np.where(ff>1,1,ff)
    return ff

