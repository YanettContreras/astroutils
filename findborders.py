import pyfits as pf
import numpy as n

#namei='G332.317-00.558_hcop_MEAN.fits'
#namef='G332.225-00.535_hcop_MEAN.fits'

def border_saboca(namei, namef):
    hdri=pf.getheader(namei)
    hdrf=pf.getheader(namef)
    nxi = hdri["NAXIS1"]
    delti = hdri["CDELT1"]
    rai = delti*(n.arange(nxi)+1-hdri["CRPIX1"]) + hdri["CRVAL1"]
    nxi2 = hdri["NAXIS2"]
    delti2 = hdri["CDELT2"]
    deci = delti2*(n.arange(nxi2)+1-hdri["CRPIX2"]) + hdri["CRVAL2"]
    nxf = hdrf["NAXIS1"]
    deltf = hdrf["CDELT1"]
    raf = deltf*(n.arange(nxf)+1-hdrf["CRPIX1"]) + hdrf["CRVAL1"]
    nxf2 = hdrf["NAXIS2"]
    deltf2 = hdrf["CDELT2"]
    decf = deltf2*(n.arange(nxf2)+1-hdrf["CRPIX2"]) + hdrf["CRVAL2"]
    ramax=n.min([raf[-1],rai[-1],raf[0],rai[0]])
    ramin=n.max([raf[-1],rai[-1],raf[0],rai[0]])
    decmax=n.min([decf[-1],deci[-1],decf[0],deci[0]])
    decmin=n.max([decf[-1],deci[-1],decf[0],deci[0]])
    crval1n=(ramax+ramin)/2.
    crval2n=(decmax+decmin)/2.
    crpix1n=n.ceil((crval1n-ramin)/deltf)
    crpix2n=n.ceil(n.abs((crval2n-decmin)/deltf))
    cdelt1n=delti
    cdelt2n=delti2
    naxis1n=n.ceil((ramax-ramin)/deltf)
    naxis2n=n.ceil(n.abs((decmax-decmin)/deltf))
    val1=[crval1n, crpix1n, cdelt1n, naxis1n]
    val2=[crval2n, crpix2n, cdelt2n, naxis2n]
    return val1,val2

def border_malt90(namei, namef):
    hdri=pf.getheader(namei)
    hdrf=pf.getheader(namef)
    nxi = hdri["NAXIS1"]
    delti = hdri["CDELT1"]
    rai = delti*(n.arange(nxi)+1-hdri["CRPIX1"]) + hdri["CRVAL1"]
    nxi2 = hdri["NAXIS2"]
    delti2 = hdri["CDELT2"]
    deci = delti2*(n.arange(nxi2)+1-hdri["CRPIX2"]) + hdri["CRVAL2"]
    nxf = hdrf["NAXIS1"]
    deltf = hdrf["CDELT1"]
    raf = deltf*(n.arange(nxf)+1-hdrf["CRPIX1"]) + hdrf["CRVAL1"]
    nxf2 = hdrf["NAXIS2"]
    deltf2 = hdrf["CDELT2"]
    decf = deltf2*(n.arange(nxf2)+1-hdrf["CRPIX2"]) + hdrf["CRVAL2"]
    naxis3 = hdri["NAXIS3"]
    cdelt3 = hdri["CDELT3"]/1000.
    crpix3 = hdri["CRPIX3"]
    crval3 = hdri["CRVAL3"]/1000.
    ramax=n.min([raf[-1],rai[-1],raf[0],rai[0]])
    ramin=n.max([raf[-1],rai[-1],raf[0],rai[0]])
    decmax=n.min([decf[-1],deci[-1],decf[0],deci[0]])
    decmin=n.max([decf[-1],deci[-1],decf[0],deci[0]])
    # print ramax,raf[-1]
    # print ramin,rai[0]
    # print decmax,decf[-1]
    # print decmin,deci[0]
    # ramax=raf[-1]
    # ramin=rai[0]
    # decmax=decf[-1]
    # decmin=deci[0]
    crval1n=(ramax+ramin)/2.
    crval2n=(decmax+decmin)/2.
    crpix1n=n.ceil((crval1n-ramin)/deltf)
    crpix2n=n.ceil(n.abs((crval2n-decmin)/deltf))
    cdelt1n=delti
    cdelt2n=delti2
    naxis1n=n.ceil((ramax-ramin)/deltf)
    naxis2n=n.ceil(n.abs((decmax-decmin)/deltf))
    val1=[crval1n, crpix1n, cdelt1n, naxis1n]
    val2=[crval2n, crpix2n, cdelt2n, naxis2n]
    val3=[crval3, crpix3, cdelt3, naxis3]
    return val1,val2,val3


