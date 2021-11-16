#include "auto_f2c.h"

inline double x(const doublereal *w, int n, int N){
  return w[(n+N)%N];
}
inline double y(const doublereal *w, int n, int N){
  return w[N+(n+N)%N];
}
inline double u(const doublereal *w, int n, int N){
  if((n+N)%N==N-1){
    return 1;
  }
  return w[2*N+(n+N)%N];
}
inline double v(const doublereal *w, int n, int N){
  if((n+N)%N==N-1){
    return 0;
  }
  return w[3*N-1+(n+N)%N];
}

int func (integer ndim, const doublereal *w, const integer *icp,
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

  for(j=0; j<N; j++){
    f[j]=-y(w,j,N)*(omega+beta*(x(w,j,N)*v(w,j,N)-u(w,j,N)*y(w,j,N)-y(w, N-1, N))+sigma*(x(w,j,N)*v(w,j+1,N)-u(w,j+1,N)*y(w,j,N)+x(w,j,N)*v(w,j-1,N)-u(w,j-1,N)*y(w,j,N)-y(w, N-2, N)-y(w, 0, N)))+gamma*(1-x(w,j,N)*x(w,j,N)-y(w,j,N)*y(w,j,N))*x(w,j,N);
    f[N+j]=x(w,j,N)*(omega+beta*(x(w,j,N)*v(w,j,N)-u(w,j,N)*y(w,j,N)-y(w, N-1, N))+sigma*(x(w,j,N)*v(w,j+1,N)-u(w,j+1,N)*y(w,j,N)+x(w,j,N)*v(w,j-1,N)-u(w,j-1,N)*y(w,j,N)-y(w, N-2, N)-y(w, 0, N)))+gamma*(1-x(w,j,N)*x(w,j,N)-y(w,j,N)*y(w,j,N))*y(w,j,N);
    if(j<N-1){
      f[2*N+j]=-v(w,j,N)*(beta*(u(w,j,N)*y(w,j,N)-x(w,j,N)*v(w,j,N)-y(w, N-1, N))+sigma*(u(w,j,N)*y(w,j+1,N)-x(w,j+1,N)*v(w,j,N)+u(w,j,N)*y(w,j-1,N)-x(w,j-1,N)*v(w,j,N)-y(w, N-2, N)-y(w, 0, N)))+gamma*(1-u(w,j,N)*u(w,j,N)-v(w,j,N)*v(w,j,N))*u(w,j,N);
      f[3*N-1+j]=u(w,j,N)*(beta*(u(w,j,N)*y(w,j,N)-x(w,j,N)*v(w,j,N)-y(w, N-1, N))+sigma*(u(w,j,N)*y(w,j+1,N)-x(w,j+1,N)*v(w,j,N)+u(w,j,N)*y(w,j-1,N)-x(w,j-1,N)*v(w,j,N)-y(w, N-2, N)-y(w, 0, N)))+gamma*(1-u(w,j,N)*u(w,j,N)-v(w,j,N)*v(w,j,N))*v(w,j,N);
    }
  }


  if (ijac == 0) {
    return 0;
  }


  // ARRAY2D(dfdu,j,j)=

  if (ijac == 1) {
    return 0;
  }

  // ARRAY2D(dfdp,j,0)=


  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par)
{
  par[0] = 0.33;
  int N = (ndim+2)/4;
  double phi0=-asin(1.0/(2*(par[0]+0.25)));
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


  doublereal order=0,order2=0,order3=0,t=0,eta=0;
  doublereal dt,weight,csum,ssum;
  int unstable=0;
  dt = getp("DTM",1,u);
  if(dt==0.0){
    eta=0;
    csum=1.0/(2*N)*(1.0+cos(eta)*u[3*(N-1)]+sin(eta)*u[3*(N-1)+N]);
    ssum=1.0/(2*N)*(0.0+cos(eta)*u[3*(N-1)+N]-sin(eta)*u[3*(N-1)]);
    for (int k=0; k<N-1; k++){
      csum+=1.0/(2*N)*(cos((k+1)*eta)*u[k]+cos(k*eta)*u[2*(N-1)+k]-sin((k+1)*eta)*u[(N-1)+k]-sin(k*eta)*u[2*(N-1)+N+k]);
      ssum+=1.0/(2*N)*(cos((k+1)*eta)*u[(N-1)+k]+cos(k*eta)*u[2*(N-1)+N+k]+sin((k+1)*eta)*u[k]+sin(k*eta)*u[2*(N-1)+k]);
    }
    order=(csum*csum+ssum*ssum);

    eta=2*3.14159265359/N;
    csum=1.0/(2*N)*(1.0+cos(eta)*u[3*(N-1)]+sin(eta)*u[3*(N-1)+N]);
    ssum=1.0/(2*N)*(0.0+cos(eta)*u[3*(N-1)+N]-sin(eta)*u[3*(N-1)]);
    for (int k=0; k<N-1; k++){
      csum+=1.0/(2*N)*(cos((k+1)*eta)*u[k]+cos(k*eta)*u[2*(N-1)+k]-sin((k+1)*eta)*u[(N-1)+k]-sin(k*eta)*u[2*(N-1)+N+k]);
      ssum+=1.0/(2*N)*(cos((k+1)*eta)*u[(N-1)+k]+cos(k*eta)*u[2*(N-1)+N+k]+sin((k+1)*eta)*u[k]+sin(k*eta)*u[2*(N-1)+k]);
    }
    order2=(csum*csum+ssum*ssum);

    eta=-2*3.14159265359/N;
    csum=1.0/(2*N)*(1.0+cos(eta)*u[3*(N-1)]+sin(eta)*u[3*(N-1)+N]);
    ssum=1.0/(2*N)*(0.0+cos(eta)*u[3*(N-1)+N]-sin(eta)*u[3*(N-1)]);
    for (int k=0; k<N-1; k++){
      csum+=1.0/(2*N)*(cos((k+1)*eta)*u[k]+cos(k*eta)*u[2*(N-1)+k]-sin((k+1)*eta)*u[(N-1)+k]-sin(k*eta)*u[2*(N-1)+N+k]);
      ssum+=1.0/(2*N)*(cos((k+1)*eta)*u[(N-1)+k]+cos(k*eta)*u[2*(N-1)+N+k]+sin((k+1)*eta)*u[k]+sin(k*eta)*u[2*(N-1)+k]);
    }
    order3=(csum*csum+ssum*ssum);
    par[10]=0;
    for(int i=1; i<ndim; i++){
      if(getp("EIG",2*i+1,u)>0){
        unstable++;
      }
    }
  }

  else{
    order=0,order2=0,order3=0;
    for (int i=0; i<NTST; i++){
      dt = getp("DTM",i+1,u);
      for (int j=0; j<NCOL+1; j++){
        weight = getp("WINT",j,u);

        eta=0;
        csum=1.0/(2*N)*(1.0+cos(eta)*ARRAY2D(u,3*(N-1),NCOL*i+j)+sin(eta)*ARRAY2D(u,3*(N-1)+N,NCOL*i+j));
        ssum=1.0/(2*N)*(0.0+cos(eta)*ARRAY2D(u,3*(N-1)+N,NCOL*i+j)-sin(eta)*ARRAY2D(u,3*(N-1),NCOL*i+j));
        for (int k=0; k<N-1; k++){
          csum+=1.0/(2*N)*(cos((k+1)*eta)*ARRAY2D(u,k,NCOL*i+j)+cos(k*eta)*ARRAY2D(u,2*(N-1)+k,NCOL*i+j)-sin((k+1)*eta)*ARRAY2D(u,(N-1)+k,NCOL*i+j)-sin(k*eta)*ARRAY2D(u,2*(N-1)+N+k,NCOL*i+j));
          ssum+=1.0/(2*N)*(cos((k+1)*eta)*ARRAY2D(u,(N-1)+k,NCOL*i+j)+cos(k*eta)*ARRAY2D(u,2*(N-1)+N+k,NCOL*i+j)+sin((k+1)*eta)*ARRAY2D(u,k,NCOL*i+j)+sin(k*eta)*ARRAY2D(u,2*(N-1)+k,NCOL*i+j));
        }
        order+=dt*weight*(csum*csum+ssum*ssum);

        eta=2*3.14159265359/N;
        csum=1.0/(2*N)*(1.0+cos(eta)*ARRAY2D(u,3*(N-1),NCOL*i+j)+sin(eta)*ARRAY2D(u,3*(N-1)+N,NCOL*i+j));
        ssum=1.0/(2*N)*(0.0+cos(eta)*ARRAY2D(u,3*(N-1)+N,NCOL*i+j)-sin(eta)*ARRAY2D(u,3*(N-1),NCOL*i+j));
        for (int k=0; k<N-1; k++){
          csum+=1.0/(2*N)*(cos((k+1)*eta)*ARRAY2D(u,k,NCOL*i+j)+cos(k*eta)*ARRAY2D(u,2*(N-1)+k,NCOL*i+j)-sin((k+1)*eta)*ARRAY2D(u,(N-1)+k,NCOL*i+j)-sin(k*eta)*ARRAY2D(u,2*(N-1)+N+k,NCOL*i+j));
          ssum+=1.0/(2*N)*(cos((k+1)*eta)*ARRAY2D(u,(N-1)+k,NCOL*i+j)+cos(k*eta)*ARRAY2D(u,2*(N-1)+N+k,NCOL*i+j)+sin((k+1)*eta)*ARRAY2D(u,k,NCOL*i+j)+sin(k*eta)*ARRAY2D(u,2*(N-1)+k,NCOL*i+j));
        }
        order2+=dt*weight*(csum*csum+ssum*ssum);

        eta=-2*3.14159265359/N;
        csum=1.0/(2*N)*(1.0+cos(eta)*ARRAY2D(u,3*(N-1),NCOL*i+j)+sin(eta)*ARRAY2D(u,3*(N-1)+N,NCOL*i+j));
        ssum=1.0/(2*N)*(0.0+cos(eta)*ARRAY2D(u,3*(N-1)+N,NCOL*i+j)-sin(eta)*ARRAY2D(u,3*(N-1),NCOL*i+j));
        for (int k=0; k<N-1; k++){
          csum+=1.0/(2*N)*(cos((k+1)*eta)*ARRAY2D(u,k,NCOL*i+j)+cos(k*eta)*ARRAY2D(u,2*(N-1)+k,NCOL*i+j)-sin((k+1)*eta)*ARRAY2D(u,(N-1)+k,NCOL*i+j)-sin(k*eta)*ARRAY2D(u,2*(N-1)+N+k,NCOL*i+j));
          ssum+=1.0/(2*N)*(cos((k+1)*eta)*ARRAY2D(u,(N-1)+k,NCOL*i+j)+cos(k*eta)*ARRAY2D(u,2*(N-1)+N+k,NCOL*i+j)+sin((k+1)*eta)*ARRAY2D(u,k,NCOL*i+j)+sin(k*eta)*ARRAY2D(u,2*(N-1)+k,NCOL*i+j));
        }
        order3+=dt*weight*(csum*csum+ssum*ssum);
      }
    }
    for(int i=1; i<ndim; i++){
      if(getp("EIG",2*i+1,u)*getp("EIG",2*i+1,u)+getp("EIG",2*i+2,u)*getp("EIG",2*i+2,u)>0.97){
        unstable++;
      }
    }
  }
  par[1]=pow(order,0.5);
  par[2]=pow(order2,0.5);
  par[3]=pow(order3,0.5);
  par[4]=ndim-unstable;
  par[5]=getp("STP",0,u);

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
