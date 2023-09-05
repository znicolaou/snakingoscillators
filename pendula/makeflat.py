#!/usr/bin/env python
import numpy as np
import argparse

if __name__ == "__main__":

	#Command line arguments
	parser = argparse.ArgumentParser(description='Driven pendula.')
	parser.add_argument("--file", type=str, required=True, dest='file', help='Base string for file output')
	parser.add_argument("--num", type=int, default=16, dest='num', help='Number of pendula')
	parser.add_argument("--dt", type=float, default=0.01, dest='dt', help='Time step')
	parser.add_argument("--invert", type=int, default=0, dest='invert', help='Zero for stable, 1 for unstable')
	args = parser.parse_args()

	dt=args.dt
	N=args.num
	n=int(1/dt)
	times=2*np.pi*dt*np.arange(0,n)[:,np.newaxis]
	xs=np.ones((n,N))
	ys=np.zeros((n,N))
	ps=np.zeros((n,N))
	us=np.ones((n,N))
	vs=np.zeros((n,N))
	qs=np.zeros((n,N))
	Xs=np.cos(times)
	Ys=np.sin(times)

	if(args.invert==1):
		xs=-xs
		us=-uso
	np.savetxt(args.file,np.concatenate([times,xs,ys,ps,us,vs,qs,Xs,Ys],axis=1))
