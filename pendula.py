#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import numpy as np
import timeit
import os
from scipy.integrate import solve_ivp

def func(t, y, N, initcycle, amp0, amp, freq0, freq, damp, spring, lengths):
	q=y[:N]
	p=y[N:]
	ampt=amp
	freqt=freq
	if t < 2*np.pi*initcycle:
		ampt=amp0+t/(2*np.pi*initcycle)*(amp-amp0)
		freqt=freq0+t/(2*np.pi*initcycle)*(freq-freq0)

	return np.concatenate( [p/lengths/freqt, (-damp*p - (1+ampt*freqt**2*np.cos(t))*np.sin(q) +spring*np.roll(lengths,1)*np.sin(np.roll(q,1)-q)+spring*np.roll(lengths,-1)*np.sin(np.roll(q,-1)-q)+spring*(np.roll(lengths,1)+np.roll(lengths,-1)-2*lengths)*np.sin(q))/freqt] )

def runsim(y0, lengths, cycles, outcycle, dt, amp, freq, initcycle=0, amp0=0, freq0=0, damp=0.1, spring=1, rtol=1e-6, atol=1e-6):
	N=len(lengths)
	sol=solve_ivp(func, [0,2*np.pi*cycles], y0, method='RK45', args=(N, initcycle, amp0, amp, freq0, freq, damp, spring, lengths), t_eval=2*np.pi*dt*np.arange(outcycle/dt,cycles/dt), rtol=rtol, atol=atol)
	ys=sol.y.T.copy()
	ys[:,:N]=np.mod(ys[:,:N]+np.pi,2*np.pi)-np.pi
	return ys

if __name__ == "__main__":

	#Command line arguments
	parser = argparse.ArgumentParser(description='Driven pendula.')
	parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for file output')
	parser.add_argument("--num", type=int, default=32, dest='num', help='Number of pendula')
	parser.add_argument("--frequency", type=float, default=3.4, dest='freq', help='Driving frequency')
	parser.add_argument("--initfrequency", type=float, default=3.4, dest='freq0', help='Initial driving frequency')
	parser.add_argument("--initamplitude", type=float, default=0.05, dest='amp0', help='Driving amplitude')
	parser.add_argument("--amplitude", type=float, default=0.05, dest='amp', help='Driving amplitude')
	parser.add_argument("--delta", type=float, default=0.5, dest='delta', help='Alternating pendulum length scale')
	parser.add_argument("--cycles", type=float, default=1000, dest='cycles', help='Simulation time in driving cycles')
	parser.add_argument("--initcycle", type=float, default=0, dest='initcycle', help='Simulation time in driving cycles')
	parser.add_argument("--outcycle", type=float, default=0, dest='outcycle', help='Cycle to start outputting')
	parser.add_argument("--dt", type=float, default=0.05, dest='dt', help='Time between outputs in driving cycles')
	parser.add_argument("--seed", type=int, default=1, dest='seed', help='Seed for random initial conditions')
	parser.add_argument("--damp", type=float, default=0.1, dest='damp', help='Damping coefficient')
	parser.add_argument("--spring", type=float, default=1.0, dest='spring', help='Spring coefficient')
	parser.add_argument("--init", type=float, default=0.01, dest='init', help='Initial random scale')
	parser.add_argument("--rtol", type=float, default=1e-6, dest='rtol', help='Relative error tolerance')
	parser.add_argument("--atol", type=float, default=1e-6, dest='atol', help='Absolute error tolerance')
	parser.add_argument("--verbose", type=int, default=1, dest='verbose', help='Verbose output')
	args = parser.parse_args()

	start = timeit.default_timer()
	N=args.num
	np.random.seed(args.seed)

	y0=np.zeros(2*N)
	if (os.path.isfile(args.filebase+"ic.npy")):
		if args.verbose==1:
			print("using initial conditions from file",flush=True)
		y0=np.load(args.filebase+"ic.npy")
	else:
		if args.verbose==1:
			print("using random initial contions",flush=True)
		y0[:N] = args.init*2*np.pi*(np.random.random(N)-0.5)
	lengths=np.array([1+args.delta*(-1)**i for i in range(N)])

	ys=runsim(y0, lengths, args.cycles, args.outcycle, args.dt, args.amp, args.freq,  initcycle=args.initcycle, amp0=args.amp0, freq0=args.freq0, damp=args.damp, spring=args.spring, rtol=args.rtol, atol=args.atol)

	#we may want to featurize the dominant frequencies instead
	norms=np.sort(np.linalg.norm(ys[:,:N],axis=0)/len(ys)**0.5)
	stop = timeit.default_timer()

	file=open(args.filebase+'out.dat','w')
	print(*sys.argv,file=file)
	print("%i %f %f %f %f %f"%(args.num, args.cycles, args.outcycle, args.dt, args.amp, args.freq), file=file)
	print(args.seed, *norms, file=file)
	print('runtime: %f' % (stop - start), file=file)
	if args.verbose==1:
		print('runtime: %f' % (stop - start))

	if args.verbose==1:
		np.save(args.filebase+"dat",ys)
		np.save(args.filebase+"fs",ys[-1])

	file.close()
