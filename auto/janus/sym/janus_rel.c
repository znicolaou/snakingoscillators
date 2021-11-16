#include "auto_f2c.h"

inline double x(const doublereal *w, int n){
  if()
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

  const doublereal *x=w;
  const doublereal *y=w+N;
  const doublereal *u=w+2*N;
  const doublereal *v=w+3*N-1;

  f[0]=-y[0]*(omega+beta*(x[0]*v[0]-u[0]*y[0]-y[N-1])+sigma*(x[0]*v[1]-u[1]*y[0]-y[j]-y[N-2]-y[0]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
  f[N+0]=x[j]*(omega+beta*(x[j]*v[j]-u[j]*y[j]-y[N-1])+sigma*(x[j]*v[j+1]-u[j+1]*y[j]+x[j]*v[j-1]-u[j-1]*y[j]-y[N-2]-y[0]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
  f[2*N+0]=-v[j]*(beta*(-x[j]*v[j]+u[j]*y[j]-y[N-1])+sigma*(u[j]*y[j+1]-x[j+1]*v[j]+u[j]*y[j-1]-x[j-1]*v[j]-y[N-2]-y[0]))+gamma*(1-u[j]*u[j]-v[j]*v[j])*u[j];
  f[3*N-1+0]=u[j]*(beta*(u[j]*y[j]-x[j]*v[j]-y[N-1])+sigma*(u[j]*y[j+1]-x[j+1]*v[j]+u[j]*y[j-1]-x[j-1]*v[j]-y[N-2]-y[0]))+gamma*(1-u[j]*u[j]-v[j]*v[j])*v[j];
  for(j=1; j<N-2; j++){
    f[j]=-y[j]*(omega+beta*(x[j]*v[j]-u[j]*y[j]-y[N-1])+sigma*(x[j]*v[j+1]-u[j+1]*y[j]+x[j]*v[j-1]-u[j-1]*y[j]-y[N-2]-y[0]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
    f[N+j]=x[j]*(omega+beta*(x[j]*v[j]-u[j]*y[j]-y[N-1])+sigma*(x[j]*v[j+1]-u[j+1]*y[j]+x[j]*v[j-1]-u[j-1]*y[j]-y[N-2]-y[0]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
    f[2*N+j]=-v[j]*(beta*(-x[j]*v[j]+u[j]*y[j]-y[N-1])+sigma*(u[j]*y[j+1]-x[j+1]*v[j]+u[j]*y[j-1]-x[j-1]*v[j]-y[N-2]-y[0]))+gamma*(1-u[j]*u[j]-v[j]*v[j])*u[j];
    f[3*N-1+j]=u[j]*(beta*(u[j]*y[j]-x[j]*v[j]-y[N-1])+sigma*(u[j]*y[j+1]-x[j+1]*v[j]+u[j]*y[j-1]-x[j-1]*v[j]-y[N-2]-y[0]))+gamma*(1-u[j]*u[j]-v[j]*v[j])*v[j];
  }


  if (ijac == 0) {
    return 0;
  }


  for(j=1; j<N-2; j++){
    ARRAY2D(dfdu,j,j)=-y[j]*(beta*(v[j+1])+sigma*(v[(j+2)]))+gamma*(1-3*x[j]*x[j]-y[j]*y[j]); //Jxx
    ARRAY2D(dfdu,j,N-1+j)=-(beta*(x[j]*v[j+1]-2*w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)]-2*w[(j+2)]*y[j]-v[1]))+gamma*(-2*y[j])*x[j]; //Jxy
    ARRAY2D(dfdu,j,2*(N-1)+j+1)=-y[j]*(beta*(-y[j])); //Jxw
    ARRAY2D(dfdu,j,2*(N-1)+(j+2)%N)=-y[j]*(sigma*(-y[j])); //Jxw
    ARRAY2D(dfdu,j,2*(N-1)+N+0)=-y[j]*(-beta); //Jxv
    ARRAY2D(dfdu,j,2*(N-1)+N+1)=-y[j]*(-sigma); //Jxv
    ARRAY2D(dfdu,j,2*(N-1)+N+j+1)=-y[j]*(beta*(x[j])); //Jxv
    ARRAY2D(dfdu,j,2*(N-1)+N+(j+2)%N)=-y[j]*(sigma*(x[j])); //Jxv

    ARRAY2D(dfdu,N-1+j,j)=(beta*(2*x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(2*x[j]*v[(j+2)]-w[(j+2)]*y[j]-v[1]))+gamma*(-2*x[j])*y[j]; //Jyx
    ARRAY2D(dfdu,N-1+j,N-1+j)=x[j]*(beta*(-w[j+1])+sigma*(-w[(j+2)]))+gamma*(1-x[j]*x[j]-3*y[j]*y[j]); //Jyy
    ARRAY2D(dfdu,N-1+j,2*(N-1)+j+1)=x[j]*(beta*(-y[j])); //Jyw
    ARRAY2D(dfdu,N-1+j,2*(N-1)+(j+2))=x[j]*(sigma*(-y[j])); //Jyw
    ARRAY2D(dfdu,N-1+j,2*(N-1)+N+0)=x[j]*(-beta); //Jyv
    ARRAY2D(dfdu,N-1+j,2*(N-1)+N+1)=x[j]*(-sigma); //Jyv
    ARRAY2D(dfdu,N-1+j,2*(N-1)+N+j+1)=x[j]*(beta*(x[j])); //Jyv
    ARRAY2D(dfdu,N-1+j,2*(N-1)+N+(j+2))=x[j]*(sigma*(x[j])); //Jyv
  }

  if (ijac == 1) {
    return 0;
  }

  for(j=0; j<N-1; j++){
    ARRAY2D(dfdp,j,0)=-y[j]*((x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]));
    ARRAY2D(dfdp,N-1+j,0)=x[j]*((x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]));
  }
  ARRAY2D(dfdp,2*(N-1)+0,0)=-v[0]*((w[0]*y[N-2]-v[0]*x[N-2]-v[1]));
  ARRAY2D(dfdp,2*(N-1)+N+0,0)=w[0]*((w[0]*y[N-2]-v[0]*x[N-2]-v[1]));
  ARRAY2D(dfdp,2*(N-1)+1,0)=-v[1]*((w[1]*0-v[1]*1-v[1]));
  ARRAY2D(dfdp,2*(N-1)+N+1,0)=w[1]*((w[1]*0-v[1]*1-v[1]));
  for(j=2; j<N; j++){
    ARRAY2D(dfdp,2*(N-1)+j,0)=-v[j]*((w[j]*y[j-2]-v[j]*x[j-2]-v[1]));
    ARRAY2D(dfdp,2*(N-1)+N+j,0)=w[j]*((w[j]*y[j-2]-v[j]*x[j-2]-v[1]));
  }

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
