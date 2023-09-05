#!/usr/bin/env python

import os
import janus
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import root
#from matplotlib import pyplot as plt

N=16
filebase='data/janus/'
num_chimeras=np.count_nonzero([st.isdigit() for st in os.listdir(filebase)])

for i in range(1,num_chimeras+1):
    os.system('mkdir -p %s/%i/steadys'%(filebase,i))
    os.system('cp {janus_rel.c,c.steady} %s/%i/steadys'%(filebase,i))
    file = open(filebase+"%i/steadys/steadys.auto"%(i), "w")
    nocycle=True
    successes=[]
    ibr=1
    while os.path.exists('%s/%i/s.fw%i_1'%(filebase,i,ibr)):
        try:
            for direction in ['fw%i'%(ibr),'rv%i'%(ibr)]:

                total=np.loadtxt(filebase+'%i/%s.dat'%(i,direction))
                if len(total)<=2:
                    continue
                if total[-1,2]>1400:
                    nocycle=False
                    print(i, ibr, direction, total[-1,0], total[-1,2])
                    omega=1.0
                    gamma=0.1
                    beta=total[-1,0]
                    sigma=total[-1,0]
                    period=total[-1,2]
                    sigma0=total[-1,0]
                    t0=0
                    sym=0

                    tf=np.load(filebase+'%i/%sfinal.dat.npy'%(i,direction))[0]
                    xf=np.load(filebase+'%i/%sfinal.dat.npy'%(i,direction))[1:N].T
                    yf=np.load(filebase+'%i/%sfinal.dat.npy'%(i,direction))[N:2*N-1].T
                    uf=np.load(filebase+'%i/%sfinal.dat.npy'%(i,direction))[2*N-1:3*N-1].T
                    vf=np.load(filebase+'%i/%sfinal.dat.npy'%(i,direction))[3*N-1:].T
                    phif=np.unwrap(np.arctan2(yf,xf),axis=0)
                    phif=np.concatenate([np.zeros(tf.shape).reshape(-1,1),phif],axis=1)
                    thetaf=np.unwrap(np.arctan2(vf,uf),axis=0)
                    Theta=1/(2*N)*(np.sum(thetaf,axis=1)+np.sum(phif,axis=1))[:,np.newaxis]
                    thetaf=np.mod(thetaf-Theta+np.pi,2*np.pi)-np.pi
                    phif=np.mod(phif-Theta+np.pi,2*np.pi)-np.pi
                    phases=np.concatenate([np.cos(phif),np.sin(phif),np.cos(thetaf),np.sin(thetaf)],axis=1)

                    norms=np.linalg.norm(np.diff(phases,axis=0)/np.diff(tf,axis=0)[:,np.newaxis],axis=1)
                    inds=find_peaks(-norms,height=-2*np.min(norms))[0]
                    successes=[]
                    stableinds=np.array(np.where(total[:,9]==4*N-2)[0])
                    stablebranches=[]
                    jumps=np.where(np.diff(stableinds)>1)[0]
                    cuts=np.zeros(2+2*len(jumps),dtype=int)
                    cuts[0]=0
                    cuts[-1]=-1
                    cuts[1:-1:2]=jumps
                    cuts[2:-1:2]=jumps+1
                    for j in range(0,len(cuts),2):
                        stablebranches=stablebranches+[stableinds[cuts[j]:cuts[j+1]]]

                    print("fw=bd([])",file=file)
                    print("rv=bd([])",file=file)

                    sols=[]
                    for ind in inds[:len(inds)//(N//2)]:
                        sol=root(lambda x: janus.janus(0,x,N,omega,sigma,beta,gamma,sigma0,0,sym),phases[ind],method='lm')
                        sols=sols+[sol]
                        successes=successes+[sol.success]
                        if sol.success:
                            z=(sol.x[:N]+1j*sol.x[N:2*N])
                            w=(sol.x[2*N:3*N]+1j*sol.x[3*N:])
                            x0=np.concatenate([np.real(z/z[0])[1:],np.imag(z/z[0])[1:],np.real(w/z[0]),np.imag(w/z[0])])
                            print("fw=fw+run('janus_rel',c='steady',IPS=1,PAR={1:%f},U={"%(sigma),end='',file=file)
                            for j in range(len(x0)):
                                print("%i:%f,"%(j+1,x0[j]),end='',file=file)
                            print('})',file=file)

                            print("rv=rv+run('janus_rel',c='steady',IPS=1,DS='-',PAR={1:%f},U={"%(sigma),end='',file=file)
                            for j in range(len(x0)):
                                print("%i:%f,"%(j+1,x0[j]),end='',file=file)
                            print('})',file=file)
                            

                    print("sv(fw,'steady%sfw')"%(direction), file=file)
                    print("sv(rv,'steady%srv')"%(direction), file=file)
                    print('total=fw+rv', file=file)
                    print("total.writeRawFilename('total%ssteady.dat')"%(direction), file=file)


        except Exception as e:
            print(e)
            pass
        ibr=ibr+1
    
    print('cl()',file=file)
    file.close()
    if not np.any(successes):
        os.system('rm '+ filebase+"%i/steadys/steadys.auto"%(i))                
