import numpy as n
import pyfits as pf
import os
import astLib as ast
import commands

def subregion(fitin=None, fittemp=None, outn=None, RADEC=False):
    """
    subregion(fitin=None, fittemp=None, outn=None)

    Crea una subregion de una imagen a partir de las coordenadas de un template
    dado, de debe ingresar un template:
    por ej:

    /malt90/mosaic/suma_G332.317-00.558/suma_G332.317-00.558_hcop.fits

    el archivo a modificar, por ej:

    /home/yanett/calan/thesis/atlasgal/mosaic/mosaico.fits

    y un nombre para el output, por ej:

    suma_G332.317-00.558_hcop
    """
    if fittemp==None:
        template = raw_input("Ingrese el archivo template : ")
    else:
        template=fittemp
    #template='/home/yanett/calan/thesis/molecular_lines/malt90/mosaic/suma_G332.317-00.558/suma_G332.317-00.558_hcop.fits'
    template2=template[:-5]+'_mom0.xy'
    print template2
    os.system('cp -r '+template2+' template2.xy')
    ftemp=pf.getdata(template)
    htemp=pf.getheader(template)

    if fitin==None:
        input1 = raw_input("Desea modificar alguno de los siguientes archivos? :\n ATLASGAL(1)\n MIPS24(2)\n GLIMPSE 3.6(3)\n GLIMPSE 4.5(4)\n GLIMPSE 5.8(5)\n GLIMPSE 8(6) OTRO(0)\n")
        if input1=='1':
            fits1= '/home/yanett/calan/thesis/atlasgal/mosaic/mosaico.fits'
        elif input1=='2':
            fits1= '/home/yanett/calan/thesis/mips/mosaico24.fits'
        elif input1=='3':
            fits1='/home/yanett/calan/thesis/glipse/mosaico3.6.fits'
        elif input1=='4':
            fits1='/home/yanett/calan/thesis/glipse/mosaico4.5.fits'
        elif input1=='5':
            fits1='/home/yanett/calan/thesis/glipse/mosaico5.8.fits'
        elif input1=='6':
            fits1='/home/yanett/calan/thesis/glipse/mosaico8.0.fits'
        else:
            fits1=raw_input("Ingrese el archivo a modificar : ")
    
    else:
        fits1=fitin
        input1='0'

    ffits=pf.getdata(fits1)
    hfits=pf.getheader(fits1)

    if outn==None:
        outname=raw_input('Ingrese el nombre del archivo de salida : ')
    else:
        outname=outn
    #outname='suma_G332.317-00.558_hcop'

    #nx1 = htemp["NAXIS1"]
    #delt1 = htemp["CDELT1"]
    #ra = delt1*(n.arange(nx1)+1-htemp["CRPIX1"]) + htemp["CRVAL1"]
    #nx2 = htemp["NAXIS2"]
    #delt2 = htemp["CDELT2"]
    #dec = delt2*(n.arange(nx2)+1-htemp["CRPIX2"]) + htemp["CRVAL2"]
    #nx3 = htemp["NAXIS3"]
    #delt3 = htemp["CDELT3"]
    #velo = (delt3*(n.arange(nx3)+1-htemp["CRPIX3"]) + htemp["CRVAL3"])/1000
    #
    #
    #nx12 = hfits["NAXIS1"]
    #delt12 = hfits["CDELT1"]
    #ra2 = delt1*(n.arange(nx12)+1-hfits["CRPIX1"]) + hfits["CRVAL1"]
    #nx22 = hfits["NAXIS2"]
    #delt22 = hfits["CDELT2"]
    #dec2 = delt2*(n.arange(nx22)+1-hfits["CRPIX2"]) + hfits["CRVAL2"]
    
    
    wcs1=ast.astWCS.WCS(template)
    wcs2=ast.astWCS.WCS(fits1)
    a=ast.astWCS.findWCSOverlap(wcs2,wcs1)
    newimag=ffits[a['wcs1Pix'][2]:a['wcs1Pix'][3],a['wcs1Pix'][1]:a['wcs1Pix'][0]]
    if RADEC==True:
        print 'USING RA-DEC MODE'
        minmax=wcs1.getImageMinMaxWCSCoords ()
        minra=ast.astCoords.convertCoords('GALACTIC','J2000',minmax[0],minmax[2],2000.)
        maxra=ast.astCoords.convertCoords('GALACTIC','J2000',minmax[1],minmax[3],2000.)
        ma1=wcs2.wcs2pix(minra[0],minra[1])
        ma2=wcs2.wcs2pix(maxra[0],maxra[1])
        a['wcs1Pix'][0]=n.min([ma1[0],ma2[0]])
        a['wcs1Pix'][1]=n.max([ma1[0],ma2[0]])
        a['wcs1Pix'][2]=n.min([ma1[1],ma2[1]])
        a['wcs1Pix'][3]=n.max([ma1[1],ma2[1]])
        print ma1,ma2,n.min(ma1[0],ma2[0]),n.max(ma1[0],ma2[0]), a
        newimag=ffits[int(a['wcs1Pix'][2]):int(a['wcs1Pix'][3]),int(a['wcs1Pix'][0]):int(a['wcs1Pix'][1])]
    
    
    centronew=wcs2.pix2wcs((a['wcs1Pix'][0]+a['wcs1Pix'][1])/2.,(a['wcs1Pix'][2]+a['wcs1Pix'][3])/2.)
    nx1new=n.shape(newimag)[1]
    nx2new=n.shape(newimag)[0]
    crpix1new=int(n.round(nx1new/2.))
    crpix2new=int(n.round(nx2new/2.))
    print nx1new,nx2new,crpix1new,crpix2new
    centronew2=wcs2.pix2wcs(int(a['wcs1Pix'][0])+crpix1new,int(a['wcs1Pix'][2])+crpix2new)
    newhea=hfits.copy()
    newhea.update("CRVAL1", centronew2[0])#crval1new)
    newhea.update("CRVAL2", centronew2[1])#crval2new)
    newhea.update("CRPIX1", crpix1new)
    newhea.update("CRPIX2", crpix2new)
    newhea.update("NAXIS1", nx1new)
    newhea.update("NAXIS2", nx2new)
    os.system('rm -r temp1.fits')
    os.system('rm -r temp1.xy')
    pf.writeto('temp1.fits', newimag, header=newhea)
    name=outname+'-'+input1
    os.system('fits in=temp1.fits op=xyin out=temp1.xy')
    os.system('rm -r '+name+".xy")
    os.system('regrid in=temp1.xy tin=template2.xy  out='+name+".xy axes=1,2" )
    os.system('fits in='+name+'.xy op=xyout out='+name+'.fits')
    #os.system('rm -r temp1.fits')
    #os.system('rm -r temp1.xy')
    #os.system('rm -r template2.xy')
    print 'Archivo creado : '+name+".xy"
    return 0

def subregion_small(dirin=None, fittemp=None, outn=None, RADEC=False,correct=True):
    """
    subregion(dirin=None, fittemp=None, outn=None)

    Crea una subregion de una imagen a partir de las coordenadas de un template
    dado, de debe ingresar un template:
    por ej:
    /malt90/mosaic/suma_G332.317-00.558/suma_G332.317-00.558_hcop.fits
    y la carpeta con los fits de donde sacar la nueva image, por ej:
    /home/yanett/calan/thesis/mips/mosaico/
    y un nombre para el output, por ej:
    suma_G332.317-00.558_hcop
    """
    ccorr=correct
    if fittemp==None:
        template = raw_input("Ingrese el archivo template : ")
    else:
        template=fittemp

    template2=template[:-5]+'_mom0.xy'
    print template2
    os.system('cp -r '+template2+' template2.xy')
 
    if dirin==None:
        dir1=raw_input("Ingrese el archivo a modificar : ")
    
    else:
        dir1=dirin
        input1='0'
        
    os.system('rm tmpfits.txt')
    os.system('ls '+dir1+'/*.fits > tmpfits.txt')
    listafits=n.loadtxt('tmpfits.txt',dtype='str')
    wcs1=ast.astWCS.WCS(template)
    bigmap=[]
    for name in listafits:
        wcs2=ast.astWCS.WCS(name)    
        minmax1=wcs1.getImageMinMaxWCSCoords()
        minmax2=wcs2.getImageMinMaxWCSCoords()
        if n.min([minmax2[0],minmax2[1]])<=minmax1[0]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[2]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
        elif n.min([minmax2[0],minmax2[1]])<=minmax1[1]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[3]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
        elif n.min([minmax2[0],minmax2[1]])<=minmax1[0]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[3]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
        elif n.min([minmax2[0],minmax2[1]])<=minmax1[1]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[2]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
    print bigmap
    print ccorr, correct
    makemosaico(bigmap, 'tempmosaico.fits', clean=True, correct=ccorr)
    if outn==None:
        outname=raw_input('Ingrese el nombre del archivo de salida : ')
    else:
        outname=outn
    name=outname+'-'+input1
    os.system('fits in=tempmosaico.fits op=xyin out=temp1.xy')
    os.system('rm -r '+name+".xy")
    os.system('regrid in=temp1.xy tin=template2.xy  out='+name+".xy axes=1,2" )
    os.system('fits in='+name+'.xy op=xyout out='+name+'.fits')
    print 'Archivo creado : '+name+".xy"
    return 0



def makemosaico(lista1, outname, clean=True, correct=True):
    print lista1
    os.system('mkdir dirtmp')
    for name in lista1:
        os.system('cp '+name+' dirtmp/.')
    os.system('mImgtbl dirtmp/ imag1.tbl')
    os.system('ls dirtmp')
    os.system('mMakeHdr imag1.tbl template1.hdr')
    os.system('mkdir projtmp')
    os.system('mProjExec -p dirtmp/ imag1.tbl template1.hdr projtmp stats1.tbl')
    os.system('mImgtbl projtmp/ images1.tbl')
    if commands.getoutput('ls -1 dirtmp | wc -l')=='1':
        print "solo un archivo"
        os.system('ls --ignore="*area*" projtmp/ > onlyone.tmp')
        onlyone=n.loadtxt('onlyone.tmp',dtype='str')
        os.system('cp projtmp/'+str(onlyone)+' '+outname)
        print 'projtmp/'+str(onlyone)
    else:
        print 'mas de dos archivos'
        os.system('mAdd -p projtmp/ images1.tbl template1.hdr uncorr.fits')
        os.system('cp uncorr.fits '+outname)
        if correct==True:
            os.system('mOverlaps images1.tbl diffs1.tbl')
            os.system('mkdir difftmp')
            os.system('mDiffExec -p projtmp/ diffs1.tbl template1.hdr difftmp')
            os.system('mFitExec diffs1.tbl fits1.tbl difftmp')
            os.system('mBgModel images1.tbl fits1.tbl correctios1.tbl')
            os.system('mkdir corrtmp')
            os.system('mBgExec -p projtmp/ images1.tbl correctios1.tbl corrtmp/')
            os.system('mAdd -p corrtmp/ images1.tbl template1.hdr mosaicotmp.fits')
            os.system('cp mosaicotmp.fits '+outname)
    if clean==True:
        os.system('rm -r projtmp')
        os.system('rm -r corrtmp')
        os.system('rm -r dirtmp')
        os.system('rm -r difftmp')
        os.system('rm uncorr.fits')
        os.system('rm mosaicotmp.fits')
        os.system('rm images1.tbl')
        os.system('rm correctios1.tbl')
        os.system('rm imag1.tbl')
        os.system('rm stats1.tbl')
        os.system('rm fits1.tbl')
        os.system('rm diffs1.tbl')
        os.system('rm template1.hdr')
    return 0

        
def subregion_small_g(dirin=None, fittemp=None, outn=None, montage=True,correct=True):
    """
    subregion(dirin=None, fittemp=None, outn=None)
    Crea una subregion de una imagen a partir de las coordenadas de un template
    dado, de debe ingresar un template:
    /home/yanett/calan/thesis/mips/mosaico/
    y un nombre para el output, por ej:
    suma_G332.317-00.558_hcop
    """
    ccorr=correct
    if fittemp==None:
        template = raw_input("Ingrese el archivo template : ")
    else:
        template=fittemp

    os.system('fits in='+fittemp+' op=xyin out=template2.xy')
    if dirin==None:
        dir1=raw_input("Ingrese el directorio de las imagenes : ")
    else:
        dir1=dirin
        input1='0'
        
    os.system('rm tmpfits.txt')
    os.system('ls '+dir1+'/*.fits > tmpfits.txt')
    listafits=n.loadtxt('tmpfits.txt',dtype='str')
    wcs1=ast.astWCS.WCS(template)
    bigmap=[]

    
    for name in listafits:
        wcs2=ast.astWCS.WCS(name)    
        minmax1=wcs1.getImageMinMaxWCSCoords()
        minmax2=wcs2.getImageMinMaxWCSCoords()
        if n.min([minmax2[0],minmax2[1]])<=minmax1[0]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[2]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
        if n.min([minmax2[0],minmax2[1]])<=minmax1[1]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[3]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
        if n.min([minmax2[0],minmax2[1]])<=minmax1[0]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[3]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
        if n.min([minmax2[0],minmax2[1]])<=minmax1[1]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[2]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
    print bigmap
    print ccorr, correct
    if montage==True:
        makemosaico(bigmap, 'tempmosaico.fits', clean=True, correct=ccorr)
    if commands.getoutput('ls -1 dirtmp | wc -l')=='1' and montage==False:
        print "Solo copiando archivo..."
        os.system('cp '+str(bigmap)+' tempmosaico.fits')
        
    if outn==None:
        outname=raw_input('Ingrese el nombre del archivo de salida : ')
    else:
        outname=outn
    name=outname+'-'+input1
    os.system('fits in=tempmosaico.fits op=xyin out=temp1.xy')
    os.system('rm -r '+name+".xy")
    os.system('regrid in=temp1.xy tin=template2.xy  out='+name+".xy axes=1,2" )
    os.system('fits in='+name+'.xy op=xyout out='+name+'.fits')
    print 'Archivo creado : '+name+".xy"
    return 0




def subregion_orig(dirin=None, fittemp=None, outn=None, montage=True,correct=True,xy=True,RADEC=False):
    """
    subregion(dirin=None, fittemp=None, outn=None)
    Crea una subregion de una imagen a partir de las coordenadas de un template
    dado, de debe ingresar un template:
    /home/yanett/calan/thesis/mips/mosaico/
    y un nombre para el output, por ej:
    suma_G332.317-00.558_hcop
    """
    ccorr=correct
    if fittemp==None:
        template = raw_input("Ingrese el archivo template : ")
    else:
        template=fittemp

    os.system('fits in='+fittemp+' op=xyin out=template2.xy')
    if dirin==None:
        dir1=raw_input("Ingrese el directorio de las imagenes : ")
    else:
        dir1=dirin
        input1='0'
        
    os.system('rm tmpfits.txt')
    os.system('ls '+dir1+'/*.fits > tmpfits.txt')
    listafits=n.loadtxt('tmpfits.txt',dtype='str')
    wcs1=ast.astWCS.WCS(template)
    bigmap=[]
    for name in listafits:
        print name
        wcs2=ast.astWCS.WCS(name)    
        minmax1=wcs1.getImageMinMaxWCSCoords()
        minmax2=wcs2.getImageMinMaxWCSCoords()
        print minmax1, minmax2
        if RADEC==True:
            minra=ast.astCoords.convertCoords('J2000','GALACTIC',minmax2[0],minmax2[2],2000.)
            maxra=ast.astCoords.convertCoords('J2000','GALACTIC',minmax2[1],minmax2[3],2000.)
            minmax2=[minra[0],maxra[0],minra[1],maxra[1]]

        if n.min([minmax2[0],minmax2[1]])<=minmax1[0]<=n.max([minmax2[0],minmax2[1]]):
            print 'true 1'
            if n.min([minmax2[2],minmax2[3]])<=minmax1[2]<=n.max([minmax2[2],minmax2[3]]):
                print 'true 2'
                bigmap=n.append(bigmap,name)
        if n.min([minmax2[0],minmax2[1]])<=minmax1[1]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[3]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
        if n.min([minmax2[0],minmax2[1]])<=minmax1[0]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[3]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
        if n.min([minmax2[0],minmax2[1]])<=minmax1[1]<=n.max([minmax2[0],minmax2[1]]):
            if n.min([minmax2[2],minmax2[3]])<=minmax1[2]<=n.max([minmax2[2],minmax2[3]]):
                bigmap=n.append(bigmap,name)
    print bigmap

    print ccorr, correct
    if montage==True:
        makemosaico(bigmap, 'tempmosaico.fits', clean=True, correct=ccorr)
    if commands.getoutput('ls -1 dirtmp | wc -l')=='1' and montage==False:
        print "Solo copiando archivo..."
        os.system('cp '+str(bigmap)+' tempmosaico.fits')
        
    if outn==None:
        outname=raw_input('Ingrese el nombre del archivo de salida : ')
    else:
        outname=outn
    name=outname+'-'+input1

    ra1=commands.getoutput('gethd in=template2.xy/crval1 format=degrees')
    dec1=commands.getoutput('gethd in=template2.xy/crval2 format=degrees')
    ra1e,dec1e=ast.astCoords.convertCoords('GALACTIC', 'J2000', float(ra1), float(dec1), 2000.)
    delta1=n.abs(float(commands.getoutput('gethd in=template2.xy/cdelt1 format=degrees')))
    xpix=float(commands.getoutput('gethd in=template2.xy/naxis1'))
    ypix=float(commands.getoutput('gethd in=template2.xy/naxis2'))
    sizey=str(delta1*xpix)
    sizex=str(delta1*ypix)
    if xy==False:
        sizey=str(delta1*ypix)
        sizex=str(delta1*xpix)
    print ra1,dec1,sizex,sizey
    #os.system('mSubimage tempmosaico.fits '+name+'.fits '+str(ra1e)+' '+str(dec1e)+' '+sizex+' '+sizey )
    os.system('mSubimage tempmosaico.fits '+name+'.fits '+str(ra1)+' '+str(dec1)+' '+sizex+' '+sizey )
    print 'Archivo creado : '+name+".fits"
    return 0




def subregion_orig2(dirin=None, fittemp=None, outn=None, montage=True,correct=True,xy=True,RADEC=False):
    """
    subregion(dirin=None, fittemp=None, outn=None)
    Crea una subregion de una imagen a partir de las coordenadas de un template
    dado, de debe ingresar un template:
    /home/yanett/calan/thesis/mips/mosaico/
    y un nombre para el output, por ej:
    suma_G332.317-00.558_hcop
    """
    ccorr=correct
    if fittemp==None:
        template = raw_input("Ingrese el archivo template : ")
    else:
        template=fittemp

    os.system('fits in='+fittemp+' op=xyin out=template2.xy')
    if dirin==None:
        dir1=raw_input("Ingrese el directorio de las imagenes : ")
    else:
        dir1=dirin
        input1='0'
        
    os.system('cp '+dir1+'/*.fits  tmp1.fits')
    name='tmp1.fits'
    wcs1=ast.astWCS.WCS(template)
    bigmap=[]
    wcs2=ast.astWCS.WCS(name)    
    minmax1=wcs1.getImageMinMaxWCSCoords()
    minmax2=wcs2.getImageMinMaxWCSCoords()
    print 'template',minmax1
    print 'archivo', minmax2
    if RADEC==True:
        minra=ast.astCoords.convertCoords('J2000','GALACTIC',minmax2[0],minmax2[2],2000.)
        maxra=ast.astCoords.convertCoords('J2000','GALACTIC',minmax2[1],minmax2[3],2000.)
        minmax2=[minra[0],maxra[0],minra[1],maxra[1]]
    print minmax1, minmax2
    
    bigmap=n.append(bigmap,name)
    print bigmap
    print ccorr, correct
    if montage==True:
        makemosaico(bigmap, 'tempmosaico.fits', clean=True, correct=ccorr)
    if commands.getoutput('ls -1 dirtmp | wc -l')=='1' and montage==False:
        print "Solo copiando archivo..."
        os.system('cp '+str(bigmap)+' tempmosaico.fits')
        
    if outn==None:
        outname=raw_input('Ingrese el nombre del archivo de salida : ')
    else:
        outname=outn
    name=outname+'-'+input1

    ra1=commands.getoutput('gethd in=template2.xy/crval1 format=degrees')
    dec1=commands.getoutput('gethd in=template2.xy/crval2 format=degrees')
    ra1e,dec1e=ast.astCoords.convertCoords('GALACTIC', 'J2000', float(ra1), float(dec1), 2000.)
    delta1=n.abs(float(commands.getoutput('gethd in=template2.xy/cdelt1 format=degrees')))
    xpix=float(commands.getoutput('gethd in=template2.xy/naxis1'))
    ypix=float(commands.getoutput('gethd in=template2.xy/naxis2'))
    sizey=str(delta1*xpix)
    sizex=str(delta1*ypix)
    if xy==False:
        sizey=str(delta1*ypix)
        sizex=str(delta1*xpix)
    print ra1e, dec1e
    os.system('mSubimage tempmosaico.fits '+name+'.fits '+str(ra1e)+' '+str(dec1e)+' '+sizex+' '+sizey )
    #os.system('mSubimage tempmosaico.fits '+name+'.fits '+str(ra1)+' '+str(dec1)+' '+sizex+' '+sizey )
    print 'Archivo creado : '+name+".fits"
    return 0


def regrig_reg(fitin=None, fittemp=None, outn=None, RADEC=False):
    """
   
    """
    if fittemp==None:
        template = raw_input("Ingrese el archivo template : ")
    else:
        template=fittemp
    
    os.system('fits in='+fittemp+' op=xyin out=template2.xy')
    ftemp=pf.getdata(fittemp)
    htemp=pf.getheader(fittemp)

    if fitin==None:
        fits1 = raw_input("Ingrese archivo input\n")
           
    else:
        fits1=fitin
        input1='0'
    os.system('fits in='+fits1+' op=xyin out=temp1.xy')
    ffits=pf.getdata(fits1)
    hfits=pf.getheader(fits1)

    if outn==None:
        outname=raw_input('Ingrese el nombre del archivo de salida : ')
    else:
        outname=outn
    
    wcs1=ast.astWCS.WCS(template)
    wcs2=ast.astWCS.WCS(fits1)
    if RADEC==True:
        print 'USING RA-DEC MODE'
        minmax=wcs1.getImageMinMaxWCSCoords ()
        minra=ast.astCoords.convertCoords('GALACTIC','J2000',minmax[0],minmax[2],2000.)
        maxra=ast.astCoords.convertCoords('GALACTIC','J2000',minmax[1],minmax[3],2000.)
    os.system('rm -r '+outname+'.xy')
    os.system('regrid in=temp1.xy tin=template2.xy  out='+outname+".xy axes=1,2" )
    os.system('fits in='+outname+'.xy op=xyout out='+outname+'.fits')
   
    print 'Archivo creado : '+outname+".fits"
    return 0
