#!/usr/bin/env python
import sys
import numpy as np
import os.path
import timeit
from scipy.sparse import csr_matrix, save_npz, load_npz, diags
from scipy.integrate import solve_ivp
from scipy.fft import fft2,ifft2
import argparse

###########################################################################
def janusrot(t, phases, N, omega, sigma, beta, gamma):
    x=np.concatenate([[1],phases[:N-1]])
    y=np.concatenate([[0],phases[N-1:2*(N-1)]])
    u=phases[2*(N-1):2*(N-1)+N]
    v=phases[2*(N-1)+N:]
    dxdt=-y*(-beta*(-x*v+y*u+v[0])-sigma*(-x*np.roll(v,1)+y*np.roll(u,1)+v[1]))+gamma*(1-x**2-y**2)*x
    dydt= x*(-beta*(-x*v+y*u+v[0])-sigma*(-x*np.roll(v,1)+y*np.roll(u,1)+v[1]))+gamma*(1-x**2-y**2)*y
    dudt=-v*(-omega-beta*(x*v-y*u+v[0])-sigma*(np.roll(x,-1)*v-np.roll(y,-1)*u+v[1]))+gamma*(1-u**2-v**2)*u
    dvdt= u*(-omega-beta*(x*v-y*u+v[0])-sigma*(np.roll(x,-1)*v-np.roll(y,-1)*u+v[1]))+gamma*(1-u**2-v**2)*v
    return np.concatenate([dxdt[1:],dydt[1:],dudt,dvdt])
###########################################################################
def runsim (N, t1, t3, dt, omega, beta, sigma, gamma, seed, filebase, output):
    start = timeit.default_timer()
    # Initial phases, from file if it exists
    np.random.seed(seed)
    phase_init = np.zeros(4*N-2,dtype=np.float64)

    if os.path.isfile(filebase + 'ic.npy'):
        phase_init = np.concatenate(np.load(filebase + 'ic.npy'))
        if output>0:
            print('using initial conditions from file',flush=True)

    elif os.path.isfile(filebase + 'ic.dat'):
        phase_init = np.fromfile(filebase + 'ic.dat',dtype=np.float64)
        if output>0:
            print('using initial conditions from file',flush=True)

    else:
        phi0=-np.pi+2*np.pi*np.random.random(N-1)
        phi1=-np.pi+2*np.pi*np.random.random(N)
        phase_init[:N-1] = np.cos(phi0)
        phase_init[N-1:2*(N-1)] = np.sin(phi0)
        phase_init[2*(N-1):2*(N-1)+N] = np.cos(phi1)
        phase_init[2*(N-1)+N:] = np.sin(phi1)
        if output>0:
            print('using random initial conditions',flush=True)

    phases=np.zeros( [int(t1/dt), 4*N-2], dtype=np.float64 )
    phase=phase_init

    sol=solve_ivp(janusrot, [0,t1], phase, method='RK45', args=(N,omega, sigma, beta, gamma), rtol=1e-6, atol=1e-6, t_eval=dt*np.arange(t1/dt))
    phases=sol.y.T.copy()
    times=sol.t

    r=np.abs(1+np.sum(phases[:,:N-1]+1j*phases[:,N-1:2*(N-1)],axis=1)+np.sum(phases[:,2*(N-1):2*(N-1)+N]+1j*phases[:,2*(N-1)+N:],axis=1))/(2*N)

    # Output
    np.save(filebase+'times.npy', times)
    np.save(filebase+'order.npy', [np.arange(0,t1,dt), r])
    np.save(filebase+'fs.npy', phases[-1, :])
    stop = timeit.default_timer()

    if output>0:
        print(filebase, np.mean(r[int(t3 / dt):]))
        print('runtime: %f' % (stop - start))

    if output>1:
        np.save(filebase+'phases.npy', phases[int(t3/dt):])

    f = open(filebase + 'out.dat', 'w')
    print(*(sys.argv), sep=' ', file=f)
    print(stop - start, np.mean(r[int(t3 / dt):]), file=f)
    f.close()
######################################################################################

if __name__ == "__main__":

    #Command line arguments
    parser = argparse.ArgumentParser(description='Numerical integration of networks of phase oscillators.')
    parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for file output.')
    parser.add_argument("--output", type=int, required=False, dest='output', default=1, help='Output style, 0 for no stdout and sparse output, 1 for stdout and sparse output, 2 for stdout and dense output. Default 1.')
    parser.add_argument("--num", type=int, required=False, dest='num', default=16, help='Number of Janus oscillators. Default 16.')
    parser.add_argument("--time", type=float, required=False, dest='time', default=10000., help='Total integration time. Detault 2000.')
    parser.add_argument("--beta", type=float, required=False, dest='beta', default=0.25, help='Internal coupling strength. Default 0.25.')
    parser.add_argument("--sigma", type=float, required=False, dest='sigma', default=0.35, help='Coupling strength. Default 1.')
    parser.add_argument("--gamma", type=float, required=False, dest='gamma', default=0.1, help='Coefficient for amplitude terms. Default 0.1.')
    parser.add_argument("--omega", type=float, required=False, dest='omega', default=1.0, help='Natural frequency gap. Default 1.0.')
    parser.add_argument("--rtime", type=float, required=False, dest='rtime', default=9000., help='Time to start averaging order parameter. Default 1000.')
    parser.add_argument("--dt", type=float, required=False, dest='dt', default=0.5, help='Time step for averaging and output. Default 0.1.')
    parser.add_argument("--seed", type=int, required=False, dest='seed', default=1, help='Initial condition random seed. Default 1.')
    args = parser.parse_args()

    N = args.num  # oscillators
    t1 = args.time  # total time
    t3 = args.rtime  # time to start averaging r
    dt = args.dt  # timestep
    beta=args.beta  # coupling strength
    sigma = args.sigma  # coupling strength
    gamma = args.gamma  # amplitude damping
    omega=args.omega #natural frequency
    seed = args.seed  # random seed
    filebase = args.filebase  # output file name
    output = args.output  # output flag


    if t3>t1:
        t3=t1-dt

    runsim(N, t1, t3, dt, omega, beta, sigma, gamma, seed, filebase, output)
