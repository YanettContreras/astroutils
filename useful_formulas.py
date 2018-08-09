import numpy as np
import pylab as plt

def mass(frequency,beta,flux,temperature,distance):
    """
        from: http://www.am.ub.edu/~robert/master/fee-part2.pdf
        """
    m1=1.6e-6
    m2=(frequency/1000)**(-(2+beta))  #in GHz
    m3=flux #in Jy
    m4=1./temperature #in K
    m5=distance**2 #in pc

    return m1*m2*m3*m4*m5


def mass_kauf(wavelegth,kappa,flux, temperature, distance):
    aa=1.439*(wavelegth**-1)*((temperature/10.0)**-1)
    m1=0.12*(np.exp(aa)-1)
    m2=(kappa/0.01)**-1
    m3=flux
    m4=(distance/100)**2
    m5=wavelegth**3
    print m1
    print m2
    return m1*m2*m3*m4*m5


def vinfall(sigma,vred,vblue,tblue,tred,tdip):
    a=sigma**2/(vred-vblue)
    b=np.log((1+np.e*(tblue/tdip))/(1+np.e*(tred/tdip)))

    return a*b

def makeplotIMF(Masses):
    # Convert to logM.
    LogMasses = np.log(np.array(Masses))

    # Plot distribution.
    plt.figure(1)
    plt.hist(Masses, 50, histtype='step', lw=3, log=True,
         range=(0.0,np.log(100.0)))
    # Overplot with Salpeter SMF.
    X = []
    Y = []
    for n in range(101):
        logM = np.log(100.0)*float(n)/100.0
        x    = np.exp(logM)
        y    = 2.0*np.power(x, 1.0-2.35)  # normalisation
        X.append(logM)
        Y.append(y)
    plt.plot(X, Y, '-', lw=3, color='black')
    plt.xlim(0.0,np.log(100.0))
    plt.xlabel(r'$\log M$', fontsize=24)
    plt.ylabel('PDF', fontsize=24)
    plt.savefig('IMF-Salpeter.png')
    plt.show()
