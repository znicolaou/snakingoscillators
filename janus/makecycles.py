#!/usr/bin/env python
import os
import numpy as np

filebase='data/randomjanus'
if not os.path.exists('%s/orders.dat'%(filebase)):
    os.system("awk '{if(FNR==4){print $0}}' %s/*out.dat > %s/orders.dat"%(filebase,filebase))
orders=np.loadtxt('%s/orders.dat'%(filebase))
orders=orders[np.argsort(orders[:,0])]

inds=np.intersect1d(np.where(orders[:,2]>10)[0],np.where(orders[:,2]<250)[0])
print(len(orders),len(inds))

b1=int(np.max(orders[inds][:,4])-np.min(orders[inds][:,4]))+1
b2=int(np.max(orders[inds][:,5])-np.min(orders[inds][:,5]))+1
binds=[1,3,4,5]
bins=[100,100,b1,b2]
counts,bins=np.histogramdd(orders[inds][:,binds],bins=bins)

chimeras=[]
periods=[]
nums=[]

for b in np.transpose(np.where(counts>=1)):
    chimerainds1=np.where(np.all([orders[inds][:,binds[i]]>=bins[i][b[i]] for i in range(4)],axis=0))[0]
    chimerainds2=np.where(np.all([orders[inds][:,binds[i]]<=bins[i][b[i]+1] for i in range(4)],axis=0))[0]
    chimerainds=np.intersect1d(chimerainds1,chimerainds2)
    chimeras=chimeras+[orders[inds[chimerainds],0].astype(int)]
    periods=periods+[int(orders[inds[chimerainds[0]],2]/0.01)]
    nums=nums+[len(chimerainds)]
print(len(chimeras))
chimerainds=np.array([np.where(orders[:,0].astype(int)==ind[0])[0][0] for ind in chimeras])
sortedchimeras=np.flip(np.lexsort(np.concatenate([orders[chimerainds][:,[1,3,2]].T,np.array(nums).reshape(1,-1)])))
last=64
sortedchimeras=sortedchimeras[:last]

seeds=''
for i in range(len(sortedchimeras)):
    seed=np.array(chimeras[sortedchimeras[i]])[0]
    seeds=seeds+str(seed)+' '
print(seeds)
print(np.transpose([np.array(nums)[sortedchimeras],np.array(periods)[sortedchimeras]]))

os.system('./sweep2.sh %s %f %f "%s"'%(filebase, 0.3, 0.3, seeds))

filebase2='data/janus'
if not os.path.exists(filebase2):
    os.mkdir(filebase2)
for i in range(len(sortedchimeras)):
    print(i+1,orders[chimerainds[sortedchimeras[i]]], nums[sortedchimeras[i]])
    phases=np.loadtxt('%s/cycles/%icycle.dat'%(filebase,np.array(chimeras[sortedchimeras[i]])[0]))

    N=int((phases.shape[1]-1+2)/4)    
    t=phases[:,0]
    x=phases[:,1:N]
    y=phases[:,N:2*N-1]
    u=phases[:,2*N-1:3*N-1]
    v=phases[:,3*N-1:]
    phi=np.unwrap(np.arctan2(y,x),axis=0)

    phi=np.concatenate([np.zeros(t.shape).reshape(-1,1),phi],axis=1)
    theta=np.unwrap(np.arctan2(v,u),axis=0)
    Theta=1/(2*N)*(np.sum(theta,axis=1)+np.sum(phi,axis=1))[:,np.newaxis]
    theta=np.mod(theta-Theta+np.pi,2*np.pi)-np.pi
    phi=np.mod(phi-Theta+np.pi,2*np.pi)-np.pi
    period=periods[sortedchimeras[i]]
    if not os.path.exists('%s/%i'%(filebase2,i+1)):
        os.mkdir('%s/%i'%(filebase2,i+1))
    os.system('cp %s/cycles/%icycle.dat %s/%i/cycle.dat'%(filebase,np.array(chimeras[sortedchimeras[i]])[0],filebase2,i+1))

