import numpy as n
import numpy
import scipy.stats as stt
import scipy.optimize as sopt

def gaussfitm(xd,yd,center,aa=8):
    lx=len(xd)
    
    data = yd[center-aa:center+aa]
    X = xd[center-aa:center+aa]
    print len(xd),len(X), len(data),center-aa, center+aa
    print center, aa, sum(abs(data))
    x = sum(X*abs(data))/sum(abs(data))
    width = n.sqrt(abs(sum((X-x)**2*abs(data))/sum(abs(data))))/2.
    max = data.max()

    fit = lambda t : max*n.exp(-(t-x)**2/(2*width**2))

    return fit, x, width, max 

def gaussfitm2(xd,yd,center,aa):
    lx=len(xd)
    
    data = yd[center-aa:center+aa]
    X = xd[center-aa:center+aa]
    print len(X), len(data)
    x = sum(X*abs(data))/sum(abs(data))
    width = n.sqrt(abs(sum((X-x)**2*abs(data))/sum(abs(data))))/2.
    maxi = data.max()

    fit = lambda t : maxi*n.exp(-(t-x)**2/(2*width**2))

    return fit, x, width, maxi



def chigauss(par,xd,yd,center,aa, cut1, cut2):
    xx,ww,mm = par
    print 'chi2',par
    ff=lambda t : mm*n.exp(-(t-xx)**2/(2*ww**2))
    if aa==0:
        aa=ww*2
    else:
        aa=aa
    # aa1=yd[:n.argmax(yd)]
    # aa2=yd[n.argmax(yd):]
    # mm=n.min([len(aa1),len(aa2)])
    # realy=(aa1[::-1][:mm]+aa2[:mm])/2

    # bb1=ff(xd)[:n.argmax(ff(xd))]
    # bb2=ff(xd)[n.argmax(ff(xd)):]
    # mm=n.min([len(bb1),len(bb2)])
    # fakef=(bb1[::-1][:mm]+bb2[:mm])/2
    # mm=n.min([len(realy),len(fakef)])
    # if mm>aa*2:
    #     mm=aa*2
    # chi=realy[:mm]-fakef[:mm]
    center1=n.argwhere(yd==n.max(yd[center-aa:center+aa]))[0][0]
    center2=n.argwhere(ff(xd)==n.max(ff(xd)[center-aa:center+aa]))[0][0]
    print center1, center2
    if center1>cut1 or center2>cut2:
        aa=5
    print n.shape(yd),n.shape(ff(xd)),center1-aa,center1+aa,center2-aa,center2+aa
    if len(yd[center1-aa:center1+aa])!=len(ff(xd)[center2-aa:center2+aa]):
        print 'algo paso'
        center2=center1
        chi=yd[center1-aa:center1+aa]-ff(xd)[center2-aa:center2+aa]
    else:
        chi=yd[center1-aa:center1+aa]-ff(xd)[center2-aa:center2+aa]

    datachi = yd[center-aa:center+aa].max()
    if n.abs(mm-datachi)>datachi*0.2:
        chi=1000
        
    print chi
    return chi

def chigaussf(par,xd,yd,center,aa, cut1, cut2):
    xx,ww = par
    mm=yd.max()
    print 'chi2',par
    ff=lambda t : mm*n.exp(-(t-xx)**2/(2*ww**2))
    if aa==0:
        aa=ww*2
    else:
        aa=aa
    center1=n.argwhere(yd==n.max(yd[center-aa:center+aa]))[0][0]
    center2=n.argwhere(ff(xd)==n.max(ff(xd)[center-aa:center+aa]))[0][0]
    print center1, center2
    if center1>cut1 or center2>cut2:
        aa=aa/2.
    print n.shape(yd),n.shape(ff(xd)),center1-aa,center1+aa,center2-aa,center2+aa
    if len(yd[center1-aa:center1+aa])!=len(ff(xd)[center2-aa:center2+aa]):
        print 'algo paso'
        center2=center1
        chi=yd[center1-aa:center1+aa]-ff(xd)[center2-aa:center2+aa]
    else:
        chi=yd[center1-aa:center1+aa]-ff(xd)[center2-aa:center2+aa]

    datachi = yd[center-aa:center+aa].max()
    if n.abs(mm-datachi)>datachi*0.2:
        chi=1000
        
    print chi
    return chi


def gaussfit(xd,yd,center, aa=10,cut1=1000, cut2=1000):
    ff, xx, ww, mm= gaussfitm (xd,yd,center, aa=aa)
    par=[xx,ww]
    try:
        res=sopt.leastsq(chigaussf, par,args=(xd,yd,center,aa,cut1,cut2), full_output=True)
        print 'lalalalallalalaa',res
        s1=res[0]
        c1=res[1]
    except TypeError:
        print 'No hay buena senal'
        s1=[0,0]
        c1=None
        mm=0
    
    print par,s1
    print 'covarianza',c1
    return [s1[0],s1[1],mm],c1

def gaussfitall(xd,yd,center, aa=10,cut1=180, cut2=80):
    ff, xx, ww, mm= gaussfitm (xd,yd,center, aa=aa)
    par=[xx,ww,mm]
    res=sopt.leastsq(chigauss, par,args=(xd,yd,center,aa,cut1,cut2), full_output=True)
    s1=res[0]
    c1=res[1]
    print par,s1
    print 'covarianza',c1
    return [s1[0],s1[1],s1[2]],c1


def gaussfitdust(xd,yd,center, aa=100,cut1=280, cut2=80):
    ff, xx, ww, mm= gaussfitm (xd,yd,center, aa=aa)
    par=[xx,ww,mm]
    s1=sopt.leastsq(chigauss, par,args=(xd,yd,center,aa,cut1,cut2))[0]
    print s1
    return [s1[0],s1[1],s1[2]]
