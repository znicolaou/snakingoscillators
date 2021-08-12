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
  int N = ndim/4;
  int j=0;

  const doublereal *x=u;
  const doublereal *y=u+N;
  const doublereal *w=u+2*N;
  const doublereal *v=u+3*N;

  for(j=0; j<N; j++){
    f[j]=-y[j]*(omega/2+beta*(x[j]*v[j]-w[j]*y[j])+sigma*(x[j]*v[(j+1)%N]-w[(j+1)%N]*y[j]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
    f[N+j]=x[j]*(omega/2+beta*(x[j]*v[j]-w[j]*y[j])+sigma*(x[j]*v[(j+1)%N]-w[(j+1)%N]*y[j]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
    f[2*N+j]=-v[j]*(-omega/2+beta*(w[j]*y[j]-v[j]*x[j])+sigma*(w[j]*y[(N+j-1)%N]-v[j]*x[(N+j-1)%N]))+gamma*(1-w[j]*w[j]-v[j]*v[j])*w[j];
    f[3*N+j]=w[j]*(-omega/2+beta*(w[j]*y[j]-v[j]*x[j])+sigma*(w[j]*y[(N+j-1)%N]-v[j]*x[(N+j-1)%N]))+gamma*(1-w[j]*w[j]-v[j]*v[j])*v[j];
  }

  if (ijac == 0) {
    return 0;
  }

  for(j=0; j<N; j++){
    ARRAY2D(dfdu,j,j)=-y[j]*(beta*v[j]+sigma*v[(j+1)%N])+gamma*(1-3*x[j]*x[j]-y[j]*y[j]); //Jxx
    ARRAY2D(dfdu,j,N+j)=-(omega/2+beta*(x[j]*v[j]-2*w[j]*y[j])+sigma*(x[j]*v[(j+1)%N]-2*w[(j+1)%N]*y[j]))+gamma*(-2*y[j])*x[j]; //Jxy
    ARRAY2D(dfdu,j,2*N+j)=-y[j]*(-beta*y[j]); //Jxw
    ARRAY2D(dfdu,j,2*N+(j+1)%N)=-y[j]*(-sigma*y[j]); //Jxw
    ARRAY2D(dfdu,j,3*N+j)=-y[j]*(beta*x[j]); //Jxv
    ARRAY2D(dfdu,j,3*N+(j+1)%N)=-y[j]*(sigma*x[j]); //Jxv

    ARRAY2D(dfdu,N+j,j)=(omega/2+beta*(2*x[j]*v[j]-w[j]*y[j])+sigma*(2*x[j]*v[(j+1)%N]-w[(j+1)%N]*y[j]))+gamma*(-2*x[j])*y[j]; //Jyx
    ARRAY2D(dfdu,N+j,N+j)=x[j]*(-beta*w[j]-sigma*w[(j+1)%N])+gamma*(1-x[j]*x[j]-3*y[j]*y[j]); //Jyy
    ARRAY2D(dfdu,N+j,2*N+j)=x[j]*(-beta*y[j]); //Jyw
    ARRAY2D(dfdu,N+j,2*N+(j+1)%N)=x[j]*(-sigma*y[j]); //Jyw
    ARRAY2D(dfdu,N+j,3*N+j)=x[j]*(beta*x[j]); //Jyv
    ARRAY2D(dfdu,N+j,3*N+(j+1)%N)=x[j]*(sigma*x[j]); //Jyv

    ARRAY2D(dfdu,2*N+j,j)=-v[j]*(-beta*v[j]); //Jwx
    ARRAY2D(dfdu,2*N+j,(N+j-1)%N)=-v[j]*(-sigma*v[j]); //Jwx
    ARRAY2D(dfdu,2*N+j,N+j)=-v[j]*(beta*w[j]); //Jwy
    ARRAY2D(dfdu,2*N+j,N+(N+j-1)%N)=-v[j]*(sigma*w[j]); //Jwy
    ARRAY2D(dfdu,2*N+j,2*N+j)=-v[j]*(beta*y[j]+sigma*y[(N+j-1)%N])+gamma*(1-3*w[j]*w[j]-v[j]*v[j]); //Jww
    ARRAY2D(dfdu,2*N+j,3*N+j)=-(-omega/2+beta*(w[j]*y[j]-2*v[j]*x[j])+sigma*(w[j]*y[(N+j-1)%N]-2*v[j]*x[(N+j-1)%N]))+gamma*(-2*v[j])*w[j]; //Jwv

    ARRAY2D(dfdu,3*N+j,j)=w[j]*(-beta*v[j]); //Jvx
    ARRAY2D(dfdu,3*N+j,(N+j-1)%N)=w[j]*(-sigma*v[j]); //Jvx
    ARRAY2D(dfdu,3*N+j,N+j)=w[j]*(beta*w[j]); //Jvy
    ARRAY2D(dfdu,3*N+j,N+(N+j-1)%N)=w[j]*(sigma*w[j]); //Jvy
    ARRAY2D(dfdu,3*N+j,2*N+j)=(-omega/2+beta*(2*w[j]*y[j]-v[j]*x[j])+sigma*(2*w[j]*y[(N+j-1)%N]-v[j]*x[(N+j-1)%N]))+gamma*(-2*w[j])*v[j]; //Jvw
    ARRAY2D(dfdu,3*N+j,3*N+j)=w[j]*(-beta*x[j]-sigma*x[(N+j-1)%N])+gamma*(1-w[j]*w[j]-3*v[j]*v[j]); //Jvv
  }

  if (ijac == 1) {
    return 0;
  }

  for(j=0; j<N; j++){
    ARRAY2D(dfdp,j,0)=-y[j]*((x[j]*v[(j+1)%N]-w[(j+1)%N]*y[j]));
    ARRAY2D(dfdp,N+j,0)=x[j]*((x[j]*v[(j+1)%N]-w[(j+1)%N]*y[j]));
    ARRAY2D(dfdp,2*N+j,0)=-v[j]*((w[j]*y[(N+j-1)%N]-v[j]*x[(N+j-1)%N]));
    ARRAY2D(dfdp,3*N+j,0)=w[j]*((w[j]*y[(N+j-1)%N]-v[j]*x[(N+j-1)%N]));
  }

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par)
{
  par[0] = 0.35;
  int N = (ndim)/4;
  double phi0=-asin(1.0/1.2);
  for (int k=0; k<N; k++){
      u[k]=1.0;
      u[N+k]=0;
      u[2*N+k]=cos(phi0);
      u[3*N+k]=sin(phi0);
  }
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
  integer u_dim1 = NDX;
  int N = NDX/4;


  doublereal order=0,t=0;
  doublereal dt,dorder,weight,csum,ssum;


  dt = getp("DTM",1,u);

  if(dt==0.0){
    csum=1.0/(2*N)*(1.0+u[3*(N-1)]);
    ssum=1.0/(2*N)*(0.0+u[3*(N-1)+N]);
    for (int k=0; k<N; k++){
      csum+=1.0/(2*N)*(u[k]+u[2*N+k]);
      ssum+=1.0/(2*N)*(u[N+k]+u[3*N+k]);
    }
    order=pow((csum*csum+ssum*ssum),0.5);
    par[10]=0;
  }

  else{
    for (int i=0; i<NTST; i++){
      dt = getp("DTM",i+1,u);
      dorder=0;
      for (int j=0; j<NCOL; j++){
        weight = getp("WINT",j,u);
        csum=0;
        ssum=0;
        for (int k=0; k<N; k++){
          csum+=1.0/(2*N)*(ARRAY2D(u,k,NCOL*i+j)+ARRAY2D(u,2*N+k,NCOL*i+j));
          ssum+=1.0/(2*N)*(ARRAY2D(u,N+k,NCOL*i+j)+ARRAY2D(u,3*N+k,NCOL*i+j));
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
