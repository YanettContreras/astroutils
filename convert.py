from Matrix import *
from RADEC import *
import math
 
def fk4tofk5(asce,decl):
    #convert coordinates from B1950 to J2000 the imput must be
    #object RA y DEC defines in RADEC.py
    #Reference: Aoki et al. AA 128,263
    
    ra = asce.RAToDegree()*0.017453292   #deg to radians
    dec = decl.DECToDegree()*0.017453292

    #if decl.sign=="-":
    #    dec=-dec
    
    A=10**(-6)*Vector([-1.62557,-0.31919,-0.13843])
    r0 = Vector([math.cos(dec)*math.cos(ra), math.cos(dec)*math.sin(ra), math.sin(dec)])

    r1 = r0-A+(r0*A)*r0

    v1= Vector([0,0,0])
    r1v1 = r1.une(v1)

    M = Matrix([[0.9999256782,0.0111820609,0.0048579479,-0.00055, 0.23849,-0.43562],[-0.0111820610,0.9999374784,-0.0000271474,-0.23854,-0.00267,0.01225],[-0.0048579477,-0.0000271765,0.9999881997,0.43574,-0.00854,0.00212],[0.0000024239502,0.0000000271066,0.0000000117766, 0.99994704, 0.01118251, 0.00485767],[-0.0000000271066,0.0000024239788,-0.0000000000658,-0.01118251,0.99995883,-0.00002714],[-0.0000000117766,-0.0000000000659,0.0000024241017,-0.00485767,-0.00002718,1.00000956]])


    rv= M * r1v1

    x= rv[0]
    y=rv[1]
    z=rv[2]


    r=math.sqrt(x**2+y**2+z**2)

    declinacion=math.asin(z/r)*57.29577951
    alpha=(math.atan2(y,x)*57.29577951)
    if alpha < 0:
        alpha = alpha+360
    else:
        alpha= alpha

    radeg=alpha
    decdeg=declinacion

    #The following part is just to convert from degrees to sexagesimal
    
    if decdeg<0:
        decsig="-"
    else:
        decsig="+"
    
    decdeg=math.fabs(decdeg)
    rahoras=int(radeg/15)
    raminut=int(((radeg/15)-rahoras)*60)
    rasegu=((((radeg/15)-rahoras)*60)-raminut)*60
        
    dechora=int(decdeg)
    decminu=int((decdeg-dechora)*60)
    decsegu=((decdeg-dechora)*60-decminu)*60
    #if int(rasegu)<10 :
     #   ra=str(rahoras)+":%2d"%raminut+":%2.3f" %rasegu
    #if int(decsegu)<10 :
     #   dec=decsig+str(dechora)+":%2d"%decminu+":%2.3f"%decsegu
    #else:
    ra=str(rahoras)+":%02d"%raminut+":%2.3f" %rasegu
    dec=decsig+str(dechora)+":%2d"%decminu+":%2.3f"%decsegu

    ascrec=RA(rahoras,raminut,rasegu)
    declin=DEC(decsig,dechora,decminu,decsegu)

    #The output is a list with coordinates in degrees(alpha, declinacion)
    #coordinates in RA DEC object (ascrec, declin) and coordinates in
    #sexagesimal formatl(ra,dec)

    return [alpha,declinacion,ascrec,declin, ra,dec]


def GALtoEC(l,b):
    #Conver galactic coordinates to Eq coordinates (Epoch 1950)
    #the imput must be latitude and longitude in degrees
    
    l=math.radians(l)
    b=math.radians(b)
    decdeg=math.asin(math.cos(b)*math.sin(l-math.radians(33))*math.sin(math.radians(62.6))+math.sin(b)*math.cos(math.radians(62.6)))
    if decdeg<0:
        decsig="-"
    else:
        decsig="+"
    A=(math.cos(b)*math.cos(l-math.radians(33))/math.cos(decdeg))
    B=(-math.sin(b)+math.sin(decdeg)*math.cos(math.radians(62.6)))/(math.sin(math.radians(62.6))*math.cos(decdeg))
    alpha1=math.radians(282.25)+math.acos(A)
    alpha2=math.radians(282.25)-math.acos(A)
    sol1=math.sin(alpha1-math.radians(282.25))
    sol2=math.sin(alpha2-math.radians(282.25))
    #print B 
    #print sol1-B
    #print sol2-B
    if math.fabs(sol1-B)<0.000000001 :

	radeg=alpha1
    if math.fabs(sol2-B)<0.000000001 :
	radeg=alpha2

    #radeg=math.radians(282.25)-math.acos(math.cos(b)*math.cos(l-math.radians(33))/math.cos(decdeg))
    decdeg=math.degrees(decdeg)
    decdeg1=math.fabs(decdeg)
    radeg=math.degrees(radeg)
    if radeg>=360:
	radeg=radeg-360
    rahoras=int(radeg/15)
    raminut=((radeg/15)-rahoras)*60
    rasegu=(raminut-int(raminut))*60    
    raminut=int(raminut)
    
    dechora=int(decdeg1)
    decminu=((decdeg1-dechora)*60)
    decsegu=(decminu-int(decminu))*60
    decminu=int(decminu)

   
    if int(rasegu)<10 :
        ra=str(rahoras)+":%02d"%raminut+":0%2.3f"%rasegu

    else:
        ra=str(rahoras)+":%02d"%raminut+":%5.3f"%rasegu
    
    if int(decsegu)<10 :
        dec=decsig+str(dechora)+":%02d"%decminu+":0%2.3f"%decsegu
   
    else:
        dec=decsig+str(dechora)+":%02d"%decminu+":%5.3f"%decsegu
    
    ascrec=RA(rahoras,raminut,rasegu)
    declin=DEC(decsig,dechora,decminu,decsegu)

    #the output is coordinates in sexagesimal format (ra,dec) and
    #coordinates in RA DEC objects (ascrec, declin)
    
    return [ra,dec,ascrec,declin,radeg,decdeg]

def ECtoGAL(ra,dec):
    #convert from B1950 equatorial coordinates to galactic coordinates
    #the imput must be an string in sexagesimal format
    
    ra=RA(float(ra[0:2]),float(ra[3:5]),float(ra[6:]))
    dec=DEC(dec[0],float(dec[1:3]),float(dec[4:6]),float(dec[7:]))

    radeg=ra.RAToDegree()
    decdeg=dec.DECToDegree()

    if dec.sign=="-":
        decdeg=-decdeg
    else:
        decdeg=decdeg
        
    rarad=math.radians(radeg)
    decrad=math.radians(decdeg)

    b=math.asin(math.sin(decrad)*math.cos(math.radians(62.6))-math.cos(decrad)*math.sin(rarad-math.radians(282.25))*math.sin(math.radians(62.6)))
    l=math.radians(33)-math.acos(math.cos(decrad)*math.cos(rarad-math.radians(282.25))/math.cos(b))

    b=math.degrees(b)
    l=math.degrees(l)

    #the output are galactic coordinates in degrees

    return [l,b]
    


    
#prueba1=RA(17,29,19.239)
#prueba2=DEC("-",36,13,25.68)
#print fk4tofk5(prueba1,prueba2)
