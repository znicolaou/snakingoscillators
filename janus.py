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

#Can we make t and phases vectors??
#This will speed up the evaluation at mesh points
###########################################################################
def janus(t, phases, N, omega, sigma, beta, gamma, sigma0, t0):
    sigmat=sigma
    if (t < t0):
        sigmat = sigma0 + (sigma - sigma0) * t / t0

    x=phases[:N]
    y=phases[N:2*N]
    u=phases[2*N:3*N]
    v=phases[3*N:4*N]
    dxdt=-y*(omega/2-beta*(-x*v+y*u)-sigmat*(-x*np.roll(v,-1)+y*np.roll(u,-1)))+gamma*(1-x**2-y**2)*x
    dydt= x*(omega/2-beta*(-x*v+y*u)-sigmat*(-x*np.roll(v,-1)+y*np.roll(u,-1)))+gamma*(1-x**2-y**2)*y
    dudt=-v*(-omega/2-beta*(x*v-y*u)-sigmat*(np.roll(x,1)*v-np.roll(y,1)*u))+gamma*(1-u**2-v**2)*u
    dvdt= u*(-omega/2-beta*(x*v-y*u)-sigmat*(np.roll(x,1)*v-np.roll(y,1)*u))+gamma*(1-u**2-v**2)*v
    return np.concatenate([dxdt,dydt,dudt,dvdt])
###########################################################################

###########################################################################
def janus_jac(t, phases, N, omega, sigma, beta, gamma, sigma0, t0):
    sigmat=sigma
    if (t < t0):
        sigmat = sigma0 + (sigma - sigma0) * t / t0

    x=phases[:N]
    y=phases[N:2*N]
    u=phases[2*N:3*N]
    v=phases[3*N:4*N]

    I1=np.identity(N)
    I2=np.roll(I1,1,axis=0)
    I3=np.roll(I1,-1,axis=0)
    # dxdt=-y*(omega/2-beta*(-x*v+y*u)-sigmat*(-x*np.roll(v,-1)+y*np.roll(u,-1)))+gamma*(1-x**2-y**2)*x
    Jxx=(-y*(-beta*(-v)-sigmat*(-np.roll(v,-1)))+gamma*(1-3*x**2-y**2))*I1
    Jxy=(-(omega/2-beta*(-x*v+2*y*u)-sigmat*(-x*np.roll(v,-1)+2*y*np.roll(u,-1)))+gamma*(-2*y)*x)*I1
    Jxu=-y*(-beta*(y*I1)-sigmat*(y*I2))
    Jxv=-y*(-beta*(-x*I1)-sigmat*(-x*I2))
    # dydt= x*(omega/2-beta*(-x*v+y*u)-sigmat*(-x*np.roll(v,-1)+y*np.roll(u,-1)))+gamma*(1-x**2-y**2)*y
    Jyx=((omega/2-beta*(-2*x*v+y*u)-sigmat*(-2*x*np.roll(v,-1)+y*np.roll(u,-1)))+gamma*(-2*x)*y)*I1
    Jyy=(x*(-beta*(u)-sigmat*(np.roll(u,-1)))+gamma*(1-x**2-3*y**2))*I1
    Jyu=x*(-beta*(y*I1)-sigmat*(y*I2))
    Jyv=x*(-beta*(-x*I1)-sigmat*(-x*I2))
    # dudt=-v*(-omega/2-beta*(x*v-y*u)-sigmat*(np.roll(x,1)*v-np.roll(y,1)*u))+gamma*(1-u**2-v**2)*u
    Jux=-v*(-beta*(v*I1)-sigmat*(v*I3))
    Juy=-v*(-beta*(-u*I1)-sigmat*(-u*I3))
    Juu=(-v*(-beta*(-y)-sigmat*(-np.roll(y,1)))+gamma*(1-3*u**2-v**2))*I1
    Juv=(-(-omega/2-beta*(2*x*v-y*u)-sigmat*(2*np.roll(x,1)*v-np.roll(y,1)*u))+gamma*(-2*v)*u)*I1
    # dvdt= u*(-omega/2-beta*(x*v-y*u)-sigmat*(np.roll(x,1)*v-np.roll(y,1)*u))+gamma*(1-u**2-v**2)*v
    Jvx=u*(-beta*(v*I1)-sigmat*(v*I3))
    Jvy=u*(-beta*(-u*I1)-sigmat*(-u*I3))
    Jvu=((-omega/2-beta*(x*v-2*y*u)-sigmat*(np.roll(x,1)*v-2*np.roll(y,1)*u))+gamma*(-2*u)*v)*I1
    Jvv=(u*(-beta*(x)-sigmat*(np.roll(x,1)))+gamma*(1-u**2-3*v**2))*I1

    return np.block([[Jxx.T,Jxy.T,Jxu.T,Jxv.T],[Jyx.T,Jyy.T,Jyu.T,Jyv.T],[Jux.T,Juy.T,Juu.T,Juv.T],[Jvx.T,Jvy.T,Jvu.T,Jvv.T]])
###########################################################################

###########################################################################
def runsim (N, t1, t3, dt, omega, beta, sigma, gamma, phase_init, sigma0=0.35, t0=-1):
    sol=solve_ivp(janus, [0,t1], phase_init, method='RK45', args=(N,omega, sigma, beta, gamma, sigma0, t0), rtol=1e-6, atol=1e-6, t_eval=dt*np.arange(t1/dt))
    phases=sol.y.T.copy()
    times=sol.t
    r=np.abs(np.sum(phases[:,:N]+1j*phases[:,N:2*N],axis=1)+np.sum(phases[:,2*N:3*N]+1j*phases[:,3*N:4*N],axis=1))/(2*N)

    return phases,times,r
######################################################################################

def order(x,y):
    return np.sum(np.diff(x)*np.abs(np.sum(y[:N]+1j*y[N:2*N],axis=0)+np.sum(y[2*N:3*N]+1j*y[3*N:],axis=0))[:-1]/(2*N))

#TODO: Calculate Floquet exponents
#TODO: Can we estimate the new limit cycle from dsigma and the Jacobian?
######################################################################################
def cont (filebase,omega,beta,gamma,sigma0,x0,y0,p0,sigmamin,sigmamax,dsigma,dsigmamax=1e-3,dsigmamin=1e-6,verbose=True, maxnodes=10000, minnodes=100, tol=1e-3, bctol=1e-3, stol=1e-4, SNum=6, coarsen=5):
    sols=[]
    sigmas=[]
    periods=[]
    orders=[]
    start=timeit.default_timer()
    N=int(len(y0)/4)
    bc=y0[0,0]
    sigma=sigma0

    def fun(ts,Xts,p):
        return p[0]*np.transpose([janus(p[0]*ts[i],Xts[:,i], N, omega, sigma, beta, gamma,sigma,-1) for i in range(len(ts))])
    def pbc(xa,xb,p):
        return np.concatenate([xb-xa,[xa[0]-bc]])
    def funjac(ts,Xts,p):
        return p[0]*np.transpose([janus_jac(p[0]*ts[i],Xts[:,i], N, omega, sigma, beta, gamma,sigma,-1) for i in range(len(ts))],(1,2,0)),  np.transpose([fun(ts,Xts,p)/p[0]],(1,0,2))
    def bcjac(xa,xb,p):
        ret0=np.zeros(4*N)
        ret1=np.zeros(4*N)
        ret0[0]=1
        return np.concatenate([-np.identity(4*N),ret0[np.newaxis,:]],axis=0),np.concatenate([np.identity(4*N),ret1[np.newaxis,:]],axis=0),np.zeros((4*N+1,1))

    start2=timeit.default_timer()
    sol=solve_bvp(fun, pbc, x0, y0, p=np.array([p0]), fun_jac=funjac, bc_jac=bcjac, max_nodes=maxnodes,tol=tol/(np.max(np.diff(x0))),bc_tol=bctol)
    stop2=timeit.default_timer()
    if verbose:
        print(sol.message,flush=True)
        print('%f\t%.3e\t%i\t%f\t%i\t%f\t'%(sigma, dsigma,len(sol.x),sol.p[0],sol.niter,stop2-start2),end='\n',flush=True)
    sols.append(sol)
    sigmas.append(sigma)
    periods.append(sol.p[0])
    orders.append(order(sol.x,sol.y))

    np.save(filebase+'_sigmas.npy',sigmas)
    np.save(filebase+'_periods.npy',periods)
    np.save(filebase+'_orders.npy',orders)
    np.save(filebase+'_times_'+str(len(sigmas)-1)+'.npy',sol.x)
    np.save(filebase+'_phases_'+str(len(sigmas)-1)+'.npy',sol.y)

    x0=sol.x
    y0=sol.y
    p0=sol.p[0]
    count=1
    SNcount=1

    while sol.success and sigma<sigmamax and sigma>sigmamin and len(x0)<=maxnodes:
        sigma=sigma+dsigma
        maxn=1.5*len(x0)
        if(sigma>sigmamax):
            sigma=sigmamax
        if(sigma<sigmamin):
            sigma=sigmamin

        try:
            start2=timeit.default_timer()
            sol=solve_bvp(fun, pbc, x0, y0, p=np.array([p0]), fun_jac=funjac, bc_jac=bcjac, max_nodes=maxn,tol=tol/(np.max(np.diff(x0))),bc_tol=bctol)
            stop2=timeit.default_timer()

            if not sol.success:
                raise Exception(sol.message)
            if (np.abs(np.max(y0)-np.max(sol.y))/np.max(y0)>5e-1):
                raise Exception('solution changed too much')
            if (len(sol.x)>maxnodes):
                raise Exception('mesh increased too much')
            if (np.abs((p0-sol.p[0])/p0)>5e-1):
                raise Exception('period changed too much '+ str(p0)+' '+str(sol.p[0]))

        except Exception as e:
            print(str(e),flush=True)
            if verbose:
                print('%f\t%.3e\t%i\t%f\t%i\t%f\t'%(sigma, dsigma,len(sol.x),sol.p[0],sol.niter,stop2-start2),end='\n',flush=True)
            sigma=sigma-dsigma
            dsigma=dsigma/2
            sol=sols[-1]
            count=1
            if np.abs(dsigma)>dsigmamin:
                continue
            else:
                if verbose:
                    print('step size too small',flush=True)
                break

        if verbose:
            print('%f\t%.3e\t%i\t%f\t%i\t%f\t'%(sigma, dsigma,len(sol.x),sol.p[0],sol.niter,stop2-start2),end='\n',flush=True)

        sols.append(sol)
        sigmas.append(sigma)
        periods.append(sol.p[0])
        orders.append(order(sol.x,sol.y))

        np.save(filebase+'_sigmas.npy',sigmas)
        np.save(filebase+'_periods.npy',periods)
        np.save(filebase+'_orders.npy',orders)
        np.save(filebase+'_times_'+str(len(sigmas)-1)+'.npy',sol.x)
        np.save(filebase+'_phases_'+str(len(sigmas)-1)+'.npy',sol.y)

        x0=sol.x
        y0=sol.y
        p0=sol.p[0]
        SNcount=SNcount+1
        count=count+1

        #Check for saddle-node
        bif=0
        if SNcount>SNum:
            ys=periods[-SNum:]
            xs=sigmas[-SNum:]
            xm=xs[-1]-(xs[-1]-xs[-3])/(ys[-1]-ys[-3])*ys[-1]
            ym=ys[-1]
            x,n=leastsq(lambda x: x[0]+x[1]*(ys-x[2])**2-xs,[xm,(xm-xs[0])/(ys[0]-ym)**2,ym])
            bif=0
            if np.abs(x[0]-sigmas[-1])<stol and (x[0]-sigmas[-1])/dsigma>0 and np.abs(x[0]-sigmas[-1])<2*np.abs(dsigma):
                bif=1
        if bif:
            count=1
            if verbose:
                print("Saddle-node expected at %f %f. Looking for second branch"%(x[0],x[2]),flush=True)
            y1=sols[-1].y
            x0=sols[-1].x
            interp=interp1d(sols[-2].x, sols[-2].y,axis=1)
            y2=np.array([interp(t) for t in x0]).T

            a=(y1*np.abs(sigmas[-2]-x[0])**0.5-y2*np.abs(sigmas[-1]-x[0])**0.5)/(np.abs(sigmas[-2]-x[0])**0.5-np.abs(sigmas[-1]-x[0])**0.5)
            b=(y2-y1)/(np.abs(sigmas[-2]-x[0])**0.5-np.abs(sigmas[-1]-x[0])**0.5)
            y0=a-b*np.abs(sigmas[-1]-x[0])**0.5
            p0=2*x[2]-periods[-1]

            start2=timeit.default_timer()
            sol2=solve_bvp(fun, pbc, x0, y0, p=np.array([p0]), fun_jac=funjac, bc_jac=bcjac, max_nodes=maxn,tol=tol/(np.max(np.diff(x0))),bc_tol=bctol)
            stop2=timeit.default_timer()


            if (not sol2.success) or (np.abs(sol2.p[0]-p0) > np.abs(sol2.p[0]-periods[-1])):
                if verbose:
                    print("Couldn't find second branch.", (not sol2.success), np.abs(sol2.p[0]-p0), np.abs(sol2.p[0]-periods[-1]), sol2.message, flush=True)
                    print('%f\t%.3e\t%i\t%f\t%i\t%f\t'%(sigma, dsigma,len(sol2.x),sol2.p[0],sol2.niter,stop2-start2),end='\n',flush=True)
                x0=sol.x
                y0=sol.y
                p0=sol.p[0]

            else:
                if np.abs(sol2.p[0]-periods[-1])/sol2.p[0] < tol/(np.max(np.diff(x0))):
                    if verbose:
                        print("Cannot distinguish branches with current tolerance",periods[-1],sol2.p[0],0,flush=True)
                    break
                if verbose:
                    print("Found second branch. Continuing.", periods[-1],sol2.p[0],flush=True)
                    print('%f\t%.3e\t%i\t%f\t%i\t%f\t'%(sigma, -dsigma,len(sol2.x),sol2.p[0],sol2.niter,stop2-start2),end='\n',flush=True)

                sols.append(sol2)
                sigmas.append(sigma)
                periods.append(sol2.p[0])
                orders.append(order(sol2.x,sol2.y))

                np.save(filebase+'_sigmas.npy',sigmas)
                np.save(filebase+'_periods.npy',periods)
                np.save(filebase+'_orders.npy',orders)
                np.save(filebase+'_times_'+str(len(sigmas)-1)+'.npy',sol2.x)
                np.save(filebase+'_phases_'+str(len(sigmas)-1)+'.npy',sol2.y)

                x0=sol2.x
                y0=sol2.y
                p0=sol2.p[0]
                dsigma=-dsigma
                SNcount=1

        if count>coarsen:
            if len(x0)>2*minnodes:
                print("Trying to coarsen.",flush=True)
                #check sol.rms_residuals to decide whether/where to coarsen?
                x1=np.concatenate([[x0[0]],x0[1:-1:2],[x0[-1]]])
                y1=np.concatenate([y0[:,:1],y0[:,1:-1:2],y0[:,-1:]],axis=1)
                start2=timeit.default_timer()
                sol2=solve_bvp(fun, pbc, x1, y1, p=np.array([p0]), fun_jac=funjac, bc_jac=bcjac, max_nodes=maxn,tol=tol/(np.max(np.diff(x1))),bc_tol=bctol)
                stop2=timeit.default_timer()
                print('%f\t%.3e\t%i\t%f\t%i\t%f\t'%(sigma, dsigma,len(sol2.x),sol2.p[0],sol2.niter,stop2-start2),end='\n',flush=True)

                if (sol2.success and len(sol2.x)<len(x0)):
                    sols.append(sol2)
                    sigmas.append(sigma)
                    periods.append(sol2.p[0])
                    orders.append(order(sol2.x,sol2.y))

                    np.save(filebase+'_sigmas.npy',sigmas)
                    np.save(filebase+'_periods.npy',periods)
                    np.save(filebase+'_orders.npy',orders)
                    np.save(filebase+'_times_'+str(len(sigmas)-1)+'.npy',sol2.x)
                    np.save(filebase+'_phases_'+str(len(sigmas)-1)+'.npy',sol2.y)

                    x0=sol2.x
                    y0=sol2.y
                    p0=sol2.p[0]

            if np.abs(dsigma)<dsigmamax:
                print("Trying to increase step.",flush=True)
                dsigma=np.sign(dsigma)*np.min([dsigmamax,np.abs(dsigma)*2])
            count=1

    return sigmas,sols
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
    parser.add_argument("--continue", type=int, required=False, dest='cont', default=0, help='Continue the solution. Default 0.')
    parser.add_argument("--dsigma", type=float, required=False, dest='dsigma', default=5e-4, help='Sigma step for continuation. Default 5e-4.')
    parser.add_argument("--dsigmamax", type=float, required=False, dest='dsigmamax', default=1e-3, help='Maximum continuation step. Default 1e-3.')
    parser.add_argument("--dsigmamin", type=float, required=False, dest='dsigmamin', default=1e-6, help='Minimum continuation step. Default 1e-6.')
    parser.add_argument("--sigmamax", type=float, required=False, dest='sigmamax', default=0.5, help='Maximum sigma for continuation. Default 0.5.')
    parser.add_argument("--sigmamin", type=float, required=False, dest='sigmamin', default=0.25, help='Minimum sigma for continuation. Default 0.25.')
    parser.add_argument("--tol", type=float, required=False, dest='tol', default=1e-4, help='Tolerance for boundary value problem. Default 1e-4.')
    parser.add_argument("--stol", type=float, required=False, dest='stol', default=1e-4, help='Tolerance for saddle-node position. Default 1e-4.')
    parser.add_argument("--coarsen", type=int, required=False, dest='coarsen', default=5, help='Number of successful steps before attempting to coarsen. Default 5.')
    parser.add_argument("--maxnodes", type=int, required=False, dest='maxnodes', default=10000, help='Maximum nodes for limit cycles. Default 10000.')
    parser.add_argument("--minnodes", type=int, required=False, dest='minnodes', default=100, help='Minimum nodes for limit cycles. Default 100.')

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
    sigmamin = args.sigmamin
    sigmamax = args.sigmamax
    dsigmamin = args.dsigmamin
    dsigmamax = args.dsigmamax
    dsigma = args.dsigma
    maxnodes = args.maxnodes
    minnodes = args.minnodes
    tol = args.tol
    stol = args.stol
    coarsen = args.coarsen

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
    phases,times,r=runsim(N, t1, t3, dt, omega, beta, sigma, gamma, phase_init)
    stop = timeit.default_timer()

    # Output
    if output>0:
        np.save(filebase+'times.npy', times)
        np.save(filebase+'order.npy', r)
        np.save(filebase+'phases.npy', phases[int(t3/dt):])
        print(filebase, np.mean(r[int(t3 / dt):]))
        print('runtime: %f' % (stop - start),flush=True)

    np.save(filebase+'fs.npy', phases[-1, :])
    f = open(filebase + 'out.dat', 'w')
    print(*(sys.argv), sep=' ', file=f)
    print(stop - start, np.mean(r[int(t3 / dt):]), file=f)
    f.close()

    if args.cont:
        phases=phases[int(t3/dt):]
        times=times[int(t3/dt):]
        minds=find_peaks(np.diff(phases[:,0]),height=0.9*np.max(np.diff(phases[:,0])))[0]
        p0=times[minds[1]]-times[minds[0]]
        x0=(times[minds[0]:minds[1]+1]-times[minds[0]])/p0
        y0=phases[minds[0]:minds[1]+1].T
        start = timeit.default_timer()
        sigmas,sols=cont(filebase+'lc_forward',omega,beta,gamma,sigma,x0,y0,p0,sigmamin,sigmamax,dsigma,dsigmamin=dsigmamin,dsigmamax=dsigmamax,maxnodes=maxnodes,minnodes=minnodes,tol=tol,bctol=tol,stol=stol,coarsen=coarsen)
        sigmas2,sols2=cont(filebase+'lc_backward',omega,beta,gamma,sigma,x0,y0,p0,sigmamin,sigmamax,-dsigma,dsigmamin=dsigmamin,dsigmamax=dsigmamax,maxnodes=maxnodes,minnodes=minnodes,tol=tol,bctol=tol,stol=stol,coarsen=coarsen)
        stop = timeit.default_timer()
        print('runtime: %f' % (stop - start),flush=True)
