import numpy as n
from column_density import *
import pyfits as pf


molecules=['n2hp','hcop','h13cop','hnc','hn13c','sio', 'hnco']

#param are mu, B , R
param=[[3.40e-18,46.586871e9,5/9.],[3.89e-18,44.594423e9,1.],[3.89e-18,43.377302e9,1.],[3.05e-18,45.331980e9,1.],[3.05e-18,43.5456e9,1.],[3.10e-18,21.711979e9,1.],[0,0,1]]

dparmol=dict(zip(molecules,param))

ratio=dict([('hcop',50.),('hnc',50.)])


def comp_coldens(mol,mom0, pos,Tex,nu,J,tau):
    ii=pf.getdata(mom0)
    print 'mu :', dparmol[mol][0]
    print 'B  :', dparmol[mol][1]
    print 'R  :', dparmol[mol][2]
    
    if mol=='sio':
        Tmb=ii[pos]*1e5  #maps are in K km/s -> K cm/s
        print Tmb
        return col_dens_sio(Tmb,dparmol[mol][1],dparmol[mol][0], Tex,J,nu)
    
    else:
        return col_dens_linmol(Tex, J, dparmol[mol][1],dparmol[mol][2],dparmol[mol][0],tau,nu)


def comp_coldens2(mol,Tex,nu,J,tau):
    print  mol
    print 'mu :', dparmol[mol][0]
    print 'B  :', dparmol[mol][1]
    print 'R  :', dparmol[mol][2]
    
    if mol=='sio':
        Tmb=tau  #maps are in K km/s -> K cm/s
        print Tmb
        return col_dens_sio(Tmb,dparmol[mol][1],dparmol[mol][0], Tex,J,nu)
    elif mol=='hnco':
        return col_dens_hnco(Tex,tau)
    
    
    else:
        return col_dens_linmol(Tex, J, dparmol[mol][1],dparmol[mol][2],dparmol[mol][0],tau,nu)




def comp_tau(mol1, mol2, file1, file2, mom1_1, mom1_2, mom2_1,mom2_2, pos, ratio, tex,j1,j2,nu1,nu2):
    
    fits1=pf.getdata(file1)[:,pos[1],pos[2]]
    fits2=pf.getdata(file2)[:,pos[1],pos[2]]
    vel1=pf.getdata(mom1_1)[pos]
    vel2=pf.getdata(mom1_2)[pos]
    width1=pf.getdata(mom2_1)[pos]
    width2=pf.getdata(mom2_2)[pos]
    
    print n.shape(fits2), vel1, width1
    
    head1=pf.getheader(file1)
    head2=pf.getheader(file2)
    
    nx3 = head1["NAXIS3"]
    delt3 = head1["CDELT3"]
    velo_1 = (delt3*(n.arange(nx3)+1-head1["CRPIX3"]) + head1["CRVAL3"])/1000
    velo1=n.argwhere(n.abs(velo_1-vel1)<0.1)[0][0]
    print 'velo and dv :',velo1, delt3
    
    nx3 = head2["NAXIS3"]
    delt3 = head2["CDELT3"]
    velo_2 = (delt3*(n.arange(nx3)+1-head2["CRPIX3"]) + head2["CRVAL3"])/1000
    velo2=n.argwhere(n.abs(velo_2-vel2)<0.1)[0][0]
    print 'velo and dv :',velo2, delt3
    
    
    if n.abs(vel1-vel2)>3.:
        print "something wrong with the velocities", vel1, vel2
    
    print delt3
    
    width=int((n.max([width1,width2])/(delt3/1000.))*3)
    print "WIDTH : ", width, velo1,velo2
    tmb1=fits1[velo1-width:velo1+width].copy()
    tmb2=fits2[velo2-width:velo2+width].copy()

    
    print len(tmb1),len(tmb2)
    taui=n.zeros([len(tmb1),2])
    for i in n.arange(len(tmb1)):
        taui[i]=solve_tau(tmb1[i], tmb2[i],ratio,tex,dparmol[mol1][1],dparmol[mol2][1],j1,j2,nu1,nu2)

    print 'Tau =', n.sum(taui)
    return fits1,fits2,tmb1,tmb2,taui





