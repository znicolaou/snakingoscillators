#include "auto_f2c.h"

int func (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
  /* System generated locals */
  integer dfdu_dim1 = ndim;
  integer dfdp_dim1 = ndim;

  double sigma = par[0];
  double omega=1.0;
  double beta=0.25;
  double gamma=1.0;
  int N = ((ndim+2)/4);
  int j=0;

  const doublereal *x=u;
  const doublereal *y=u+(N-1);
  const doublereal *w=u+2*(N-1);
  const doublereal *v=u+2*(N-1)+N;

  for(j=0; j<N-2; j++){
    f[j]=-y[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[j+2]-w[j+2]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
    f[N-1+j]=x[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[j+2]-w[j+2]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
  }

  f[N-2]=-y[N-2]*(beta*(x[N-2]*v[N-1]-w[N-1]*y[N-2]-v[0])+sigma*(x[N-2]*v[0]-w[0]*y[N-2]-v[1]))+gamma*(1-x[N-2]*x[N-2]-y[N-2]*y[N-2])*x[N-2];
  f[N-1+N-2]=x[N-2]*(beta*(x[N-2]*v[N-1]-w[N-1]*y[N-2]-v[0])+sigma*(x[N-2]*v[0]-w[0]*y[N-2]-v[1]))+gamma*(1-x[N-2]*x[N-2]-y[N-2]*y[N-2])*y[N-2];

  f[2*(N-1)+0]=-v[0]*(-omega+beta*(w[0]*0-v[0]*1-v[0])+sigma*(w[0]*y[N-2]-v[0]*x[N-2]-v[1]))+gamma*(1-w[0]*w[0]-v[0]*v[0])*w[0];
  f[2*(N-1)+N+0]=w[0]*(-omega+beta*(w[0]*0-v[0]*1-v[0])+sigma*(w[0]*y[N-2]-v[0]*x[N-2]-v[1]))+gamma*(1-w[0]*w[0]-v[0]*v[0])*v[0];

  f[2*(N-1)+1]=-v[1]*(-omega+beta*(w[1]*y[0]-v[1]*x[0]-v[0])+sigma*(w[1]*0-v[1]*1-v[1]))+gamma*(1-w[1]*w[1]-v[1]*v[1])*w[1];
  f[2*(N-1)+N+1]=w[1]*(-omega+beta*(w[1]*y[0]-v[1]*x[0]-v[0])+sigma*(w[1]*0-v[1]*1-v[1]))+gamma*(1-w[1]*w[1]-v[1]*v[1])*v[1];

  for(j=2; j<N; j++){
    f[2*(N-1)+j]=-v[j]*(-omega+beta*(w[j]*y[j-1]-v[j]*x[j-1]-v[0])+sigma*(w[j]*y[j-2]-v[j]*x[j-2]-v[1]))+gamma*(1-w[j]*w[j]-v[j]*v[j])*w[j];
    f[2*(N-1)+N+j]=w[j]*(-omega+beta*(w[j]*y[j-1]-v[j]*x[j-1]-v[0])+sigma*(w[j]*y[j-2]-v[j]*x[j-2]-v[1]))+gamma*(1-w[j]*w[j]-v[j]*v[j])*v[j];
  }

  if (ijac == 0) {
    return 0;
  }
  //Replace I1, I2, and I3 with loops over nonzero indices i=j, i=j+1, i=j-1...
  // I1=np.identity(N)
  // I2=np.roll(I1,1,axis=0)
  // I3=np.roll(I1,-1,axis=0)
  // # dxdt=-y*(omega/2-beta*(-x*v+y*u)-sigmat*(-x*np.roll(v,-1)+y*np.roll(u,-1)))+gamma*(1-x**2-y**2)*x
  // Jxx=(-y*(-beta*(-v)-sigmat*(-np.roll(v,-1)))+gamma*(1-3*x**2-y**2))*I1
  // Jxy=(-(omega/2-beta*(-x*v+2*y*u)-sigmat*(-x*np.roll(v,-1)+2*y*np.roll(u,-1)))+gamma*(-2*y)*x)*I1
  // Jxu=-y*(-beta*(y*I1)-sigmat*(y*I2))
  // Jxv=-y*(-beta*(-x*I1)-sigmat*(-x*I2))
  // # dydt= x*(omega/2-beta*(-x*v+y*u)-sigmat*(-x*np.roll(v,-1)+y*np.roll(u,-1)))+gamma*(1-x**2-y**2)*y
  // Jyx=((omega/2-beta*(-2*x*v+y*u)-sigmat*(-2*x*np.roll(v,-1)+y*np.roll(u,-1)))+gamma*(-2*x)*y)*I1
  // Jyy=(x*(-beta*(u)-sigmat*(np.roll(u,-1)))+gamma*(1-x**2-3*y**2))*I1
  // Jyu=x*(-beta*(y*I1)-sigmat*(y*I2))
  // Jyv=x*(-beta*(-x*I1)-sigmat*(-x*I2))
  // # dudt=-v*(-omega/2-beta*(x*v-y*u)-sigmat*(np.roll(x,1)*v-np.roll(y,1)*u))+gamma*(1-u**2-v**2)*u
  // Jux=-v*(-beta*(v*I1)-sigmat*(v*I3))
  // Juy=-v*(-beta*(-u*I1)-sigmat*(-u*I3))
  // Juu=(-v*(-beta*(-y)-sigmat*(-np.roll(y,1)))+gamma*(1-3*u**2-v**2))*I1
  // Juv=(-(-omega/2-beta*(2*x*v-y*u)-sigmat*(2*np.roll(x,1)*v-np.roll(y,1)*u))+gamma*(-2*v)*u)*I1
  // # dvdt= u*(-omega/2-beta*(x*v-y*u)-sigmat*(np.roll(x,1)*v-np.roll(y,1)*u))+gamma*(1-u**2-v**2)*v
  // Jvx=u*(-beta*(v*I1)-sigmat*(v*I3))
  // Jvy=u*(-beta*(-u*I1)-sigmat*(-u*I3))
  // Jvu=((-omega/2-beta*(x*v-2*y*u)-sigmat*(np.roll(x,1)*v-2*np.roll(y,1)*u))+gamma*(-2*u)*v)*I1
  // Jvv=(u*(-beta*(x)-sigmat*(np.roll(x,1)))+gamma*(1-u**2-3*v**2))*I1

  //
  // for(j=0; j<N; j++){
  //   ARRAY2D(dfdu,j,j) = 0;
  // }

  if (ijac == 1) {
    return 0;
  }

  // for(j=0; j<N; j++){
  //   ARRAY2D(dfdp,0,0) = 0;
  // }

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par)
{
  par[0] = 0.35;
  int N = (ndim+2)/4;
  double phi0=-asin(1.0/1.2);
  for (int k=0; k<N-1; k++){
      u[k]=1.0;
      u[(N-1)+k]=0;
      u[2*(N-1)+k]=cos(phi0);
      u[2*(N-1)+N+k]=sin(phi0);
  }
  u[3*(N-1)]=cos(phi0);
  u[3*(N-1)+N]=sin(phi0);
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int pvls (integer ndim, const doublereal *u,
          doublereal *par)
{
  extern double getp();
  integer NDX  = getp("NDX", 0, u);
  integer NTST  = getp("NTST", 0, u);
  integer NCOL  = getp("NCOL", 0, u);
  integer IPS  = getp("IPS", 0, u);
  integer u_dim1 = NDX;
  int N = (ndim+2)/4;


  doublereal order=0,t=0;
  doublereal dt,dorder,weight,csum,ssum;
  dt = getp("DTM",1,u);

  if(dt==0.0){
    csum=1.0/(2*N)*(1.0+u[3*(N-1)]);
    ssum=1.0/(2*N)*(0.0+u[3*(N-1)+N]);
    for (int k=0; k<N-1; k++){
      csum+=1.0/(2*N)*(u[k]+u[2*(N-1)]);
      ssum+=1.0/(2*N)*(u[(N-1)+k]+u[2*(N-1)+N+k]);
    }
    order=pow((csum*csum+ssum*ssum),0.5);
  }

  else{
    for (int i=0; i<NTST; i++){
      dt = getp("DTM",i+1,u);
      dorder=0;
      for (int j=0; j<NCOL; j++){
        weight = getp("WINT",j,u);
        csum=1.0/(2*N)*(1.0+ARRAY2D(u,3*(N-1),NCOL*i+j));
        ssum=1.0/(2*N)*(0.0+ARRAY2D(u,3*(N-1)+N,NCOL*i+j));
        for (int k=0; k<N-1; k++){
          csum+=1.0/(2*N)*(ARRAY2D(u,k,NCOL*i+j)+ARRAY2D(u,2*(N-1)+k,NCOL*i+j));
          ssum+=1.0/(2*N)*(ARRAY2D(u,(N-1)+k,NCOL*i+j)+ARRAY2D(u,2*(N-1)+N+k,NCOL*i+j));
        }
        dorder+=weight*pow((csum*csum+ssum*ssum),0.5);
      }
      t+=dt;
      order+=dt*dorder;
    }
  }
  par[1]=order;
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc)
{
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint)
{
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
