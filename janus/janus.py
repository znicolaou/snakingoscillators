#!/usr/bin/env python
import sys
import numpy as np
import os.path
import timeit
from scipy.sparse import csr_matrix, save_npz, load_npz, diags
from scipy.integrate import solve_ivp, solve_bvp
from scipy.fft import fft2,ifft2
from scipy.optimize import leastsq
from scipy.signal import argrelmax,find_peaks,argrelmin
from scipy.interpolate import interp1d
import argparse

###########################################################################
def janus(t, phases, N, omega, sigma, beta, gamma, sigma0, t0, sym):
    sigmat=sigma
    if (t < t0):
        sigmat = sigma0 + (sigma - sigma0) * t / t0

    x=phases[:N]
    y=phases[N:2*N]
    u=phases[2*N:3*N]
    v=phases[3*N:4*N]
    dxdt=-y*(omega/2-beta*(-x*v+y*u)-sigmat*(-x*np.roll(v,-1)+y*np.roll(u,-1))-sym*sigmat*(-x*np.roll(v,1)+y*np.roll(u,1)))+gamma*(1-x**2-y**2)*x
    dydt= x*(omega/2-beta*(-x*v+y*u)-sigmat*(-x*np.roll(v,-1)+y*np.roll(u,-1))-sym*sigmat*(-x*np.roll(v,1)+y*np.roll(u,1)))+gamma*(1-x**2-y**2)*y
    dudt=-v*(-omega/2-beta*(x*v-y*u)-sigmat*(np.roll(x,1)*v-np.roll(y,1)*u)-sym*sigmat*(np.roll(x,-1)*v-np.roll(y,-1)*u))+gamma*(1-u**2-v**2)*u
    dvdt= u*(-omega/2-beta*(x*v-y*u)-sigmat*(np.roll(x,1)*v-np.roll(y,1)*u)-sym*sigmat*(np.roll(x,-1)*v-np.roll(y,-1)*u))+gamma*(1-u**2-v**2)*v
    return np.concatenate([dxdt,dydt,dudt,dvdt])
###########################################################################

###########################################################################
def runsim (N, t1, t3, dt, omega, beta, sigma, gamma, sym, phase_init, sigma0=0.35, t0=-1):
    sol=solve_ivp(janus, [0,t1], phase_init, method='RK45', args=(N,omega, sigma, beta, gamma, sigma0, t0, sym), rtol=1e-6, atol=1e-6, t_eval=dt*np.arange(t1/dt))
    phases=sol.y.T.copy()
    times=sol.t

    phi=np.arctan2(phases[int(t3/dt):,N:2*N],phases[int(t3/dt):,:N])
    theta=np.arctan2(phases[int(t3/dt):,3*N:],phases[int(t3/dt):,2*N:3*N])
    angles=np.concatenate([np.mod(theta-theta[:,0][:,np.newaxis]+np.pi,2*np.pi)-np.pi,np.mod(phi-theta[:,0][:,np.newaxis]+np.pi,2*np.pi)-np.pi],axis=1)
    norms=np.linalg.norm(np.mod(angles-angles[-1]+np.pi,2*np.pi)-np.pi,axis=1)
    mins=np.array(argrelmin(norms))[0]
    mins=mins[np.where(norms[mins]<1E-1)[0]]

    n0=10
    p0=0.0
    if(len(mins)>5):
        periods=np.diff(times[mins])
        periods=periods[np.where(periods>0)[0]]
        if np.min(periods)>10.0 and np.std(periods) < 1.0:
            p0=np.min(periods)
            n0=10*int(p0/dt)

    R=np.sum(np.exp(1j*(phi[-n0:]-phi[-n0:,[0]]))+np.exp(1j*(theta[-n0:]-phi[-n0:,[0]])),axis=1)/(2*N)
    Rp=np.sum(np.exp(1j*(phi[-n0:]-phi[-n0:,[0]]+np.arange(N)/N*2*np.pi))+np.exp(1j*(theta[-n0:]-phi[-n0:,[0]]+np.arange(N)/N*2*np.pi)),axis=1)/(2*N)
    Rn=np.sum(np.exp(1j*(phi[-n0:]-phi[-n0:,[0]]-np.arange(N)/N*2*np.pi))+np.exp(1j*(theta[-n0:]-phi[-n0:,[0]]-np.arange(N)/N*2*np.pi)),axis=1)/(2*N)
    w=int(np.round(np.sum(np.diff(np.unwrap(np.angle(R))))/(10*2*np.pi)))
    wp=int(np.round(np.sum(np.diff(np.unwrap(np.angle(Rp))))/(10*2*np.pi)))
    wn=int(np.round(np.sum(np.diff(np.unwrap(np.angle(Rn))))/(10*2*np.pi)))
    r=np.mean(np.abs(R)**2,axis=0)**0.5
    rp=np.mean(np.abs(Rp)**2,axis=0)**0.5
    rn=np.mean(np.abs(Rn)**2,axis=0)**0.5

    return phases, times, r, p0, rp-rn, w, wp-wn

######################################################################################
if __name__ == "__main__":

    #Command line arguments
    parser = argparse.ArgumentParser(description='Numerical integration of networks of phase oscillators.')
    parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for file output.')
    parser.add_argument("--output", type=        # orderp=np.mean(rp[-n0:])**0.5
int, required=False, dest='output', default=1, help='Output style, 0 for no stdout and sparse output, 1 for stdout and sparse output, 2 for stdout and dense output. Default 1.')
    parser.add_argument("--num", type=int, required=False, dest='num', default=16, help='Number of Janus oscillators. Default 16.')
    parser.add_argument("--time", type=float, required=False, dest='time', default=1000., help='Total integration time. Detault 2000.')
    parser.add_argument("--beta", type=float, required=False, dest='beta', default=0.25, help='Internal coupling strength. Default 0.25.')
    parser.add_argument("--sigma", type=float, required=False, dest='sigma', default=0.35, help='Coupling strength. Default 1.')
    parser.add_argument("--gamma", type=float, required=False, dest='gamma', default=0.1, help='Coefficient for amplitude terms. Default 0.1.')
    parser.add_argument("--omega", type=float, required=False, dest='omega', default=1.0, help='Natural frequency gap. Default 1.0.')
    parser.add_argument("--rtime", type=float, required=False, dest='rtime', default=0., help='Time to start averaging order parameter. Default 1000.')
    parser.add_argument("--dt", type=float, required=False, dest='dt', default=0.1, help='Time step for averaging and output. Default 0.1.')
    parser.add_argument("--seed", type=int, required=False, dest='seed', default=1, help='Initial condition random seed. Default 1.')
    parser.add_argument("--symmetrize", type=int, required=False, dest='sym', default=1, help='Symmetrize the coupling. 1 for symmetric, 0 for chiral. Default 0.')

    args = parser.parse_args()

    N = args.num
    t1 = args.time
    t3 = args.rtime
    dt = args.dt
    beta=args.beta
    sigma = args.sigma
    gamma = args.gamma
    omega=args.omega
    seed = args.seed
    filebase = args.filebase
    output = args.output
    sym = args.sym

    if t3>t1:
        t3=t1-dt


    # Initial phases, from file if it exists
    np.random.seed(seed)
    phase_init = np.zeros(4*N,dtype=np.float64)

    if os.path.isfile(filebase + 'ic.npy'):
        phase_init = np.load(filebase + 'ic.npy')
        if output>0:
            print('using initial conditions from file',flush=True)

    elif os.path.isfile(filebase + 'ic.dat'):
        phase_init = np.fromfile(filebase + 'ic.dat',dtype=np.float64)
        if output>0:
            print('using initial conditions from file',flush=True)

    else:
        phi0=-np.pi+2*np.pi*np.random.random(N)
        phi1=-np.pi+2*np.pi*np.random.random(N)
        phase_init[:N] = np.cos(phi0)
        phase_init[N:2*N] = np.sin(phi0)
        phase_init[2*N:3*N] = np.cos(phi1)
        phase_init[3*N:] = np.sin(phi1)
        if output>0:
            print('using random initial conditions',flush=True)

    start = timeit.default_timer()
    phases, times, r, p0, r2, w, w2=runsim(N, t1, t3, dt, omega, beta, sigma, gamma, sym, phase_init)
    stop = timeit.default_timer()

    # Output
    if output>0:
        np.save(filebase+'times.npy', times)
        np.save(filebase+'phases.npy', phases[int(t3/dt):])
        print(filebase, r)
        print('runtime: %f' % (stop - start),flush=True)

    np.save(filebase+'fs.npy', phases[-1, :])
    f = open(filebase + 'out.dat', 'w')

    if(output>0):
        print(seed, r, p0, r2, w, w2)

    if(output>2):
        if p0>0:
            period=int(p0/dt)
            phi=np.arctan2(phases[-period:,N:2*N],phases[-period:,:N])
            theta=np.arctan2(phases[-period:,3*N:],phases[-period:,2*N:3*N])
            cycle=np.concatenate([np.arange(period).reshape(-1,1)*args.dt,np.cos(phi[:,1:]-phi[:,[0]]),np.sin(phi[:,1:]-phi[:,[0]]),np.cos(theta-phi[:,[0]]),np.sin(theta-phi[:,[0]])],axis=1)
            np.savetxt(filebase+'cycle.dat',cycle)

    print(*(sys.argv), sep=' ', file=f)
    print(stop - start, file=f)
    print("%i %i %f %f %f"%(N, seed, args.time, args.rtime, args.dt), sep=' ', file=f)
    print(seed,r,p0,r2,w,w2, file=f)
    f.close()
