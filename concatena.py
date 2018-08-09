import numpy as n
import scipy.io as sio
import string

def concatena(lista,out):
    lista2=n.loadtxt(lista,dtype='str')
    cat=n.empty([0,n.shape(n.loadtxt(lista2[0], dtype='str'))[1]])
    
    for file in lista2:
        print file
        cat2=n.loadtxt(file,dtype='str')
        cat=n.append(cat,cat2,axis=0)
        
    f=open(out,'w')

    print >> f, '\n'.join(map(lambda x:'\t'.join(x),cat))

    f.close()
        
   
def concatena_id(lista,out):
    lista2=n.loadtxt(lista,dtype='str')
    cat=n.empty([0,n.shape(n.loadtxt(lista2[0], dtype='str'))[1]+1])
    
    for file in lista2:
        print file
        cat2=n.loadtxt(file,dtype='str')
        cat3=n.empty([n.shape(cat2)[0],n.shape(cat2)[1]+1],dtype='S20')
        cat3[:,:n.shape(cat2)[1]]=cat2.copy()
        namefile=string.join(file.split('.')[1:3],'.')
        print namefile
        cat3[:,n.shape(cat2)[1]]=n.repeat(namefile,len(cat2))
        cat=n.append(cat,cat3,axis=0)
        
    f=open(out,'w')

    print >> f, '\n'.join(map(lambda x:'\t'.join(x),cat))

    f.close()
        
   

    
