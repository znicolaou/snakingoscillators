#include "auto_f2c.h"

int func (integer ndim, const doublereal *w, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
  /* System generated locals */
  integer dfdu_dim1 = ndim;
  integer dfdp_dim1 = ndim;

  double amp = par[0];
  double omega = par[1];
  double delta=0.35;
  double eta=0.1;
  double gamma=1;
  int N = ((ndim-2)/6);
  int j=0;

  const doublereal *x=w;
  const doublereal *y=w+N;
  const doublereal *p=w+2*N;
  const doublereal *u=w+3*N;
  const doublereal *v=w+4*N;
  const doublereal *q=w+5*N;
  const doublereal X=w[6*N];
  const doublereal Y=w[6*N+1];

  for(j=0; j<N; j++){
    f[0*N+j]=-p[j]*y[j]/(1+delta)+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
    f[1*N+j]= p[j]*x[j]/(1+delta)+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
    f[2*N+j]=-eta*p[j]-(1-omega*omega*X-4*delta)*y[j]+(1-delta)*(x[j]*v[j]-u[j]*y[j]+x[j]*v[(j-1+N)%N]-y[j]*u[(j-1+N)%N]);
    f[3*N+j]=-q[j]*v[j]/(1-delta)+gamma*(amp*amp-u[j]*u[j]-v[j]*v[j])*u[j];
    f[4*N+j]= q[j]*u[j]/(1-delta)+gamma*(amp*amp-u[j]*u[j]-v[j]*v[j])*v[j];
    f[5*N+j]=-eta*q[j]-(1-omega*omega*X+4*delta)*v[j]+(1+delta)*(u[j]*y[j]-x[j]*v[j]+u[j]*y[(j+1)%N]-v[j]*x[(j+1)%N]);
  }
  f[6*N+0]=-omega*Y+gamma*(1-X*X-Y*Y)*X;
  f[6*N+1]= omega*X+gamma*(1-X*X-Y*Y)*Y;

  if (ijac == 0) {
    return 0;
  }

  if (ijac == 1) {
    return 0;
  }

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par)
{
  par[0] = 0.1;
  par[1]=3.4;
  int j, N = (ndim-2)/6;
  for (j=0; k<N; k++){
      u[0*N+j]=1.0;
      u[1*N+j]=0;
      u[2*N+j]=0;
      u[3*N+j]=1.0;
      u[4*N+j]=0;
      u[5*N+j]=0;
  }
  u[6*N]=1;
  u[6*N+1]=0;
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
  int N = (ndim-2)/6;


  doublereal t=0;
  doublereal dt,weight;
  int unstable=0;
  dt = getp("DTM",1,u);

  if(dt==0.0){
    par[10]=0;
    for(int i=1; i<ndim; i++){
      if(getp("EIG",2*i+1,u)>0){
        unstable++;
      }
    }
  }

  else{

    for (int i=0; i<NTST; i++){
      dt = getp("DTM",i+1,u);
      for (int j=0; j<NCOL+1; j++){
        weight = getp("WINT",j,u);
      }
    }
    for(int i=1; i<ndim; i++){
      if(getp("EIG",2*i+1,u)*getp("EIG",2*i+1,u)+getp("EIG",2*i+2,u)*getp("EIG",2*i+2,u)>0.97){
        unstable++;
      }
    }
  }

  par[2]=ndim-unstable;
  par[3]=getp("STP",0,u);

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
