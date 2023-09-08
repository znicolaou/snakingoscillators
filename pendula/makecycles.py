#!/usr/bin/env python

from scipy.cluster import vq
import numpy as np
import os

if not os.path.exists('data/randompendula/orders.dat'):
    os.system("awk '{if(FNR==3){print $0}}' data/randompendula/*out.dat > data/randompendula/orders.dat")
vals=np.loadtxt('data/randompendula/orders.dat')

num=20
np.random.seed(1)
clusters=vq.kmeans(vals[:,1:],num)[0]
clusters=clusters[np.argsort(np.linalg.norm(clusters,axis=1))]
ids=vq.vq(vals[:,1:],clusters)[0]
print(np.max(ids))


chimeras=np.array([vals[np.where(ids==i)[0][0],0] for i in range(np.max(ids))],dtype=int)
print(chimeras)
print(np.unique(ids,return_counts=True)[1])

seeds=''
for i in range(len(chimeras)):
    seed=chimeras[i]
    seeds=seeds+str(seed)+' '
print(seeds)
os.system('./sweep2.sh %s'%(seeds))


if not os.path.exists('data/pendula'):
    os.mkdir('data/pendula/')
for i in range(len(chimeras)):
    ys=np.load('data/randompendula/%idat.npy'%(chimeras[i]))
    file=open('data/randompendula/%iout.dat'%(chimeras[i]))

    N,tmax,t1,dt,amp,omega=np.array(file.readlines()[1].split(),dtype=np.float64)
    N=int(N)
    delta=0.25
    omega=3.5
    phases=ys[:201,:N]
    X=np.cos(2*np.pi/100*(np.arange(201))).reshape(201,1)
    Y=np.sin(2*np.pi/100*(np.arange(201))).reshape(201,1)
    t=(np.arange(201)/200*4*np.pi).reshape(201,1)
    x=np.cos(phases[:,::2])
    y=np.sin(phases[:,::2])
    p=omega*(1+delta)*(np.roll(phases[:-1,::2],-1,axis=0)-np.roll(phases[:-1,::2],1,axis=0))/(2*2*np.pi*dt)
    p=np.concatenate([p,[p[0]]])

    u=np.cos(phases[:,1::2])
    v=np.sin(phases[:,1::2])
    q=omega*(1-delta)*(np.roll(phases[:-1,1::2],-1,axis=0)-np.roll(phases[:-1,1::2],1,axis=0))/(2*2*np.pi*dt)
    q=np.concatenate([q,[q[0]]])
    cycle=np.concatenate([t,x,y,p,u,v,q,X,Y],axis=1)
    if not os.path.exists('data/pendula/%i'%(i)):
        os.mkdir('data/pendula/%i'%(i))
    np.savetxt('data/pendula/%i/cycle.dat'%(i),cycle)
