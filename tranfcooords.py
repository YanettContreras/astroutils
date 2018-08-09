import astLib as ast
import numpy as n

def trf_coord(inputSystem, outputSystem, coordX, coordY, epoch, dec=True):
    out=n.zeros([len(coordX),2])
    print len(coordX)
    for i in range(len(coordX)):
        if dec==False:
            cdX=ast.astCoords.hms2decimal(coordX[i],":")
            cdY=ast.astCoords.dms2decimal(coordY[i],":")
            print cdX,cdY
        else:
            cdX=coordX[i]
            cdY=coordY[i]
        out[i]= ast.astCoords.convertCoords(inputSystem, outputSystem, cdX, cdY, epoch)
        

    return out


def sextodec(coordX, coordY):
    out=n.zeros([len(coordX),2])
    print len(coordX)
    for i in range(len(coordX)):
        cdX=ast.astCoords.hms2decimal(coordX[i],":")
        cdY=ast.astCoords.dms2decimal(coordY[i],":")
        out[i]=cdX,cdY
    return out
