import pyfits as pf
import numpy as n
import findborders as fbor
import os
import string as st
import numpy as n

archives=[]
arch1='T'
while(arch1!=''):
    arch1=str(raw_input('Ingrese archivo :'))
    if arch1!='':
        archives=n.append(archives, arch1)

val1,val2=fbor.border2d(name1+'.fits', name2+'.fits')
i=0
os.system('rm -r out_*.xy')
for archv in archives:
    #os.system('fits in='+archv+'.fits op=xyin out='+archv+'.xy')
    namei=archv+'.xy'
    os.system('regrid in='+namei+" out=out_"+str(i)+".xy axes=1,2 desc="+str(val1[0])+','+str(val1[1])+','+str(val1[2])+','+str(val1[3])+','+str(val2[0])+','+str(val2[1])+','+str(val2[2])+','+str(val2[3]))
    i=i+1
    
os.system('ls -d out_*.xy')
outcomb='mosaico_0.xy'

for ind in n.arange(i-1)+1:
    os.system('prthd in=out_'+str(ind)+'.xy')
    os.system('imstat in=out_'+str(ind)+'.xy |tail -2\n')
    rms0=n.append(rms0, str(raw_input(chr(27)+"[0;33m"+'Valor de rms '+str(ind)+': '+chr(27)+"[0m")))
    outcomb=n.append(outcomb,'out_'+str(ind)+'.xy')

rms1=st.join(rms0,',')
output=st.join(outcomb,',')
os.system('rm -r suma.xy')
os.system('imcomb in='+output+' rms='+rms1+' out=suma.xy')
os.system('imsub in=suma.xy out=suma_'+archives[0]+'_.xy')
os.system('addhis in=suma_'+archives[0]+'_.xy comment='+st.join(archives,','))


