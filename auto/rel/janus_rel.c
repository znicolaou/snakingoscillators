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

  for(j=0; j<N-1; j++){
    f[j]=-y[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
    f[N-1+j]=x[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
  }
  f[2*(N-1)+0]=-v[0]*(-omega+beta*(-2*v[0])+sigma*(w[0]*y[N-2]-v[0]*x[N-2]-v[1]))+gamma*(1-w[0]*w[0]-v[0]*v[0])*w[0];
  f[2*(N-1)+N+0]=w[0]*(-omega+beta*(-2*v[0])+sigma*(w[0]*y[N-2]-v[0]*x[N-2]-v[1]))+gamma*(1-w[0]*w[0]-v[0]*v[0])*v[0];
  f[2*(N-1)+1]=-v[1]*(-omega+beta*(w[1]*y[0]-v[1]*x[0]-v[0])+sigma*(-2*v[1]))+gamma*(1-w[1]*w[1]-v[1]*v[1])*w[1];
  f[2*(N-1)+N+1]=w[1]*(-omega+beta*(w[1]*y[0]-v[1]*x[0]-v[0])+sigma*(-2*v[1]))+gamma*(1-w[1]*w[1]-v[1]*v[1])*v[1];
  for(j=2; j<N; j++){
    f[2*(N-1)+j]=-v[j]*(-omega+beta*(w[j]*y[j-1]-v[j]*x[j-1]-v[0])+sigma*(w[j]*y[j-2]-v[j]*x[j-2]-v[1]))+gamma*(1-w[j]*w[j]-v[j]*v[j])*w[j];
    f[2*(N-1)+N+j]=w[j]*(-omega+beta*(w[j]*y[j-1]-v[j]*x[j-1]-v[0])+sigma*(w[j]*y[j-2]-v[j]*x[j-2]-v[1]))+gamma*(1-w[j]*w[j]-v[j]*v[j])*v[j];
  }

  if (ijac == 0) {
    return 0;
  }

  //just make loop j=1; j<N-2 and add two cases outside loop...
  //j=0
  // f[j]=-y[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
  ARRAY2D(dfdu,0,0)=-y[0]*(beta*(v[1])+sigma*(v[2]))+gamma*(1-3*x[0]*x[0]-y[0]*y[0]); //Jxx
  ARRAY2D(dfdu,0,N-1+0)=-(beta*(x[0]*v[1]-2*w[1]*y[0]-v[0])+sigma*(x[0]*v[2]-2*w[2]*y[0]-v[1]))+gamma*(-2*y[0])*x[0]; //Jxy
  ARRAY2D(dfdu,0,2*(N-1)+0+1)=-y[0]*(beta*(-y[0])); //Jxw
  ARRAY2D(dfdu,0,2*(N-1)+0+2)=-y[0]*(sigma*(-y[0])); //Jxw
  ARRAY2D(dfdu,0,2*(N-1)+N+0)=-y[0]*(-beta); //Jxv
  ARRAY2D(dfdu,0,2*(N-1)+N+1)=-y[0]*(beta*(x[0])-sigma); //Jxv
  ARRAY2D(dfdu,0,2*(N-1)+N+0+2)=-y[0]*(sigma*(x[0])); //Jxv

  // f[N-1+j]=x[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
  ARRAY2D(dfdu,N-1+0,0)=(beta*(2*x[0]*v[1]-w[1]*y[0]-v[0])+sigma*(2*x[0]*v[2]-w[2]*y[0]-v[1]))+gamma*(-2*x[0])*y[0]; //Jyx
  ARRAY2D(dfdu,N-1+0,N-1+0)=x[0]*(beta*(-w[1])+sigma*(-w[2]))+gamma*(1-x[0]*x[0]-3*y[0]*y[0]); //Jyy
  ARRAY2D(dfdu,N-1+0,2*(N-1)+0+1)=x[0]*(beta*(-y[0])); //Jyw
  ARRAY2D(dfdu,N-1+0,2*(N-1)+0+2)=x[0]*(sigma*(-y[0])); //Jyw
  ARRAY2D(dfdu,N-1+0,2*(N-1)+N+0)=x[0]*(-beta); //Jyv
  ARRAY2D(dfdu,N-1+0,2*(N-1)+N+1)=x[0]*(beta*(x[0])-sigma); //Jyv
  ARRAY2D(dfdu,N-1+0,2*(N-1)+N+0+2)=x[0]*(sigma*(x[0])); //Jyv

  //j=N-2
  // f[j]=-y[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
  ARRAY2D(dfdu,N-2,N-2)=-y[N-2]*(beta*(v[N-1])+sigma*(v[0]))+gamma*(1-3*x[N-2]*x[N-2]-y[N-2]*y[N-2]); //Jxx
  ARRAY2D(dfdu,N-2,N-1+N-2)=-(beta*(x[N-2]*v[N-1]-2*w[N-1]*y[N-2]-v[0])+sigma*(x[N-2]*v[0]-2*w[0]*y[N-2]-v[1]))+gamma*(-2*y[N-2])*x[N-2]; //Jxy
  ARRAY2D(dfdu,N-2,2*(N-1)+N-2+1)=-y[N-2]*(beta*(-y[N-2])); //Jxw
  ARRAY2D(dfdu,N-2,2*(N-1)+0)=-y[N-2]*(sigma*(-y[N-2])); //Jxw
  ARRAY2D(dfdu,N-2,2*(N-1)+N+0)=-y[N-2]*(-beta+sigma*(x[N-2])); //Jxv
  ARRAY2D(dfdu,N-2,2*(N-1)+N+1)=-y[N-2]*(-sigma); //Jxv
  ARRAY2D(dfdu,N-2,2*(N-1)+N+N-1)=-y[N-2]*(beta*(x[N-2])); //Jxv

  // f[N-1+j]=x[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
  ARRAY2D(dfdu,N-1+N-2,N-2)=(beta*(2*x[N-2]*v[N-1]-w[N-1]*y[N-2]-v[0])+sigma*(2*x[N-2]*v[0]-w[0]*y[N-2]-v[1]))+gamma*(-2*x[N-2])*y[N-2]; //Jyx
  ARRAY2D(dfdu,N-1+N-2,N-1+N-2)=x[N-2]*(beta*(-w[N-1])+sigma*(-w[0]))+gamma*(1-x[N-2]*x[N-2]-3*y[N-2]*y[N-2]); //Jyy
  ARRAY2D(dfdu,N-1+N-2,2*(N-1)+N-1)=x[N-2]*(beta*(-y[N-2])); //Jyw
  ARRAY2D(dfdu,N-1+N-2,2*(N-1)+0)=x[N-2]*(sigma*(-y[N-2])); //Jyw
  ARRAY2D(dfdu,N-1+N-2,2*(N-1)+N+0)=x[N-2]*(-beta+sigma*(x[N-2])); //Jyv
  ARRAY2D(dfdu,N-1+N-2,2*(N-1)+N+1)=x[N-2]*(-sigma); //Jyv
  ARRAY2D(dfdu,N-1+N-2,2*(N-1)+N+N-1)=x[N-2]*(beta*(x[N-2])); //Jyv

  for(j=1; j<N-2; j++){
    // f[j]=-y[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
    ARRAY2D(dfdu,j,j)=-y[j]*(beta*(v[j+1])+sigma*(v[(j+2)]))+gamma*(1-3*x[j]*x[j]-y[j]*y[j]); //Jxx
    ARRAY2D(dfdu,j,N-1+j)=-(beta*(x[j]*v[j+1]-2*w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)]-2*w[(j+2)]*y[j]-v[1]))+gamma*(-2*y[j])*x[j]; //Jxy
    ARRAY2D(dfdu,j,2*(N-1)+j+1)=-y[j]*(beta*(-y[j])); //Jxw
    ARRAY2D(dfdu,j,2*(N-1)+(j+2)%N)=-y[j]*(sigma*(-y[j])); //Jxw
    ARRAY2D(dfdu,j,2*(N-1)+N+0)=-y[j]*(-beta); //Jxv
    ARRAY2D(dfdu,j,2*(N-1)+N+1)=-y[j]*(-sigma); //Jxv
    ARRAY2D(dfdu,j,2*(N-1)+N+j+1)=-y[j]*(beta*(x[j])); //Jxv
    ARRAY2D(dfdu,j,2*(N-1)+N+(j+2)%N)=-y[j]*(sigma*(x[j])); //Jxv

    // f[N-1+j]=x[j]*(beta*(x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(x[j]*v[(j+2)%N]-w[(j+2)%N]*y[j]-v[1]))+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
    ARRAY2D(dfdu,N-1+j,j)=(beta*(2*x[j]*v[j+1]-w[j+1]*y[j]-v[0])+sigma*(2*x[j]*v[(j+2)]-w[(j+2)]*y[j]-v[1]))+gamma*(-2*x[j])*y[j]; //Jyx
    ARRAY2D(dfdu,N-1+j,N-1+j)=x[j]*(beta*(-w[j+1])+sigma*(-w[(j+2)]))+gamma*(1-x[j]*x[j]-3*y[j]*y[j]); //Jyy
    ARRAY2D(dfdu,N-1+j,2*(N-1)+j+1)=x[j]*(beta*(-y[j])); //Jyw
    ARRAY2D(dfdu,N-1+j,2*(N-1)+(j+2))=x[j]*(sigma*(-y[j])); //Jyw
    ARRAY2D(dfdu,N-1+j,2*(N-1)+N+0)=x[j]*(-beta); //Jyv
    ARRAY2D(dfdu,N-1+j,2*(N-1)+N+1)=x[j]*(-sigma); //Jyv
    ARRAY2D(dfdu,N-1+j,2*(N-1)+N+j+1)=x[j]*(beta*(x[j])); //Jyv
    ARRAY2D(dfdu,N-1+j,2*(N-1)+N+(j+2))=x[j]*(sigma*(x[j])); //Jyv
  }
  //j=0
  // f[2*(N-1)+0]=-v[0]*(-omega+beta*(-2*v[0])+sigma*(w[0]*y[N-2]-v[0]*x[N-2]-v[1]))+gamma*(1-w[0]*w[0]-v[0]*v[0])*w[0];
  ARRAY2D(dfdu,2*(N-1)+0,N-2)=-v[0]*(sigma*(-v[0])); //Jwx
  ARRAY2D(dfdu,2*(N-1)+0,N-1+N-2)=-v[0]*(sigma*(w[0])); //Jwy
  ARRAY2D(dfdu,2*(N-1)+0,2*(N-1)+0)=-v[0]*(sigma*(y[N-2]))+gamma*(1-3*w[0]*w[0]-v[0]*v[0]); //Jww
  ARRAY2D(dfdu,2*(N-1)+0,2*(N-1)+N+0)=-(-omega+beta*(-4*v[0])+sigma*(w[0]*y[N-2]-2*v[0]*x[N-2]-v[1]))+gamma*(-2*v[0])*w[0]; //Jwv
  ARRAY2D(dfdu,2*(N-1)+0,2*(N-1)+N+1)=-v[0]*(-sigma); //Jwv
  // f[2*(N-1)+N+0]=w[0]*(-omega+beta*(-2*v[0])+sigma*(w[0]*y[N-2]-v[0]*x[N-2]-v[1]))+gamma*(1-w[0]*w[0]-v[0]*v[0])*v[0];
  ARRAY2D(dfdu,2*(N-1)+N+0,N-2)=w[0]*(sigma*(-v[0])); //Jvx
  ARRAY2D(dfdu,2*(N-1)+N+0,N-1+N-2)=w[0]*(sigma*(w[0])); //Jvy
  ARRAY2D(dfdu,2*(N-1)+N+0,2*(N-1)+0)=(-omega+beta*(-2*v[0])+sigma*(2*w[0]*y[N-2]-v[0]*x[N-2]-v[1]))+gamma*(-2*w[0])*v[0]; //Jvw
  ARRAY2D(dfdu,2*(N-1)+N+0,2*(N-1)+N+0)=w[0]*(beta*(-2)+sigma*(-x[N-2]))+gamma*(1-w[0]*w[0]-3*v[0]*v[0]); //Jvv
  ARRAY2D(dfdu,2*(N-1)+N+0,2*(N-1)+N+1)=w[0]*(-sigma); //Jvv

  //j=1
  // f[2*(N-1)+1]=-v[1]*(-omega+beta*(w[1]*y[0]-v[1]*x[0]-v[0])+sigma*(-2*v[1]))+gamma*(1-w[1]*w[1]-v[1]*v[1])*w[1];
  ARRAY2D(dfdu,2*(N-1)+1,0)=-v[1]*(beta*(-v[1])); //Jwx
  ARRAY2D(dfdu,2*(N-1)+1,N-1+0)=-v[1]*(beta*(w[1])); //Jwy
  ARRAY2D(dfdu,2*(N-1)+1,2*(N-1)+1)=-v[1]*(beta*(y[0]))+gamma*(1-3*w[1]*w[1]-v[1]*v[1]); //Jww
  ARRAY2D(dfdu,2*(N-1)+1,2*(N-1)+N+1)=-(-omega+beta*(w[1]*y[0]-2*v[1]*x[0]-v[0])+sigma*(-4*v[1]))+gamma*(-2*v[1])*w[1]; //Jwv
  ARRAY2D(dfdu,2*(N-1)+1,2*(N-1)+N+0)=-v[1]*(-beta); //Jwv
  // f[2*(N-1)+N+1]=w[1]*(-omega+beta*(w[1]*y[0]-v[1]*x[0]-v[0])+sigma*(-2*v[1]))+gamma*(1-w[1]*w[1]-v[1]*v[1])*v[1];
  ARRAY2D(dfdu,2*(N-1)+N+1,0)=w[1]*(beta*(-v[1])); //Jvx
  ARRAY2D(dfdu,2*(N-1)+N+1,N-1+0)=w[1]*(beta*(w[1])); //Jvy
  ARRAY2D(dfdu,2*(N-1)+N+1,2*(N-1)+1)=(-omega+beta*(2*w[1]*y[0]-v[1]*x[0]-v[0])+sigma*(-2*v[1]))+gamma*(-2*w[1])*v[1]; //Jvw
  ARRAY2D(dfdu,2*(N-1)+N+1,2*(N-1)+N+1)=w[1]*(beta*(-x[0])+sigma*(-2))+gamma*(1-w[1]*w[1]-3*v[1]*v[1]); //Jvv
  ARRAY2D(dfdu,2*(N-1)+N+1,2*(N-1)+N+0)=w[1]*(-beta); //Jvv

  for(j=2; j<N; j++){
    // f[2*(N-1)+j]=-v[j]*(-omega+beta*(w[j]*y[j-1]-v[j]*x[j-1]-v[0])+sigma*(w[j]*y[j-2]-v[j]*x[j-2]-v[1]))+gamma*(1-w[j]*w[j]-v[j]*v[j])*w[j];
    ARRAY2D(dfdu,2*(N-1)+j,j-1)=-v[j]*(beta*(-v[j])); //Jwx
    ARRAY2D(dfdu,2*(N-1)+j,j-2)=-v[j]*(sigma*(-v[j])); //Jwx
    ARRAY2D(dfdu,2*(N-1)+j,N-1+j-1)=-v[j]*(beta*(w[j])); //Jwy
    ARRAY2D(dfdu,2*(N-1)+j,N-1+j-2)=-v[j]*(sigma*(w[j])); //Jwy
    ARRAY2D(dfdu,2*(N-1)+j,2*(N-1)+j)=-v[j]*(beta*(y[j-1])+sigma*(y[j-2]))+gamma*(1-3*w[j]*w[j]-v[j]*v[j]); //Jww
    ARRAY2D(dfdu,2*(N-1)+j,2*(N-1)+N+j)=-(-omega+beta*(w[j]*y[j-1]-2*v[j]*x[j-1]-v[0])+sigma*(w[j]*y[j-2]-2*v[j]*x[j-2]-v[1]))+gamma*(-2*v[j])*w[j]; //Jwv
    ARRAY2D(dfdu,2*(N-1)+j,2*(N-1)+N+0)=-v[j]*(-beta); //Jwv
    ARRAY2D(dfdu,2*(N-1)+j,2*(N-1)+N+1)=-v[j]*(-sigma); //Jwv
    // f[2*(N-1)+N+j]=w[j]*(-omega+beta*(w[j]*y[j-1]-v[j]*x[j-1]-v[0])+sigma*(w[j]*y[j-2]-v[j]*x[j-2]-v[1]))+gamma*(1-w[j]*w[j]-v[j]*v[j])*v[j];
    ARRAY2D(dfdu,2*(N-1)+N+j,j-1)=w[j]*(beta*(-v[j])); //Jvx
    ARRAY2D(dfdu,2*(N-1)+N+j,j-2)=w[j]*(sigma*(-v[j])); //Jvx
    ARRAY2D(dfdu,2*(N-1)+N+j,N-1+j-1)=w[j]*(beta*(w[j])); //Jvy
    ARRAY2D(dfdu,2*(N-1)+N+j,N-1+j-2)=w[j]*(sigma*(w[j])); //Jvy
    ARRAY2D(dfdu,2*(N-1)+N+j,2*(N-1)+j)=(-omega+beta*(2*w[j]*y[j-1]-v[j]*x[j-1]-v[0])+sigma*(2*w[j]*y[j-2]-v[j]*x[j-2]-v[1]))+gamma*(-2*w[j])*v[j]; //Jvw
    ARRAY2D(dfdu,2*(N-1)+N+j,2*(N-1)+N+j)=w[j]*(beta*(-x[j-1])+sigma*(-x[j-2]))+gamma*(1-w[j]*w[j]-3*v[j]*v[j]); //Jvv
    ARRAY2D(dfdu,2*(N-1)+N+j,2*(N-1)+N+0)=w[j]*(-beta); //Jvv
    ARRAY2D(dfdu,2*(N-1)+N+j,2*(N-1)+N+1)=w[j]*(-sigma); //Jvv
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


  doublereal order=0,t=0;
  doublereal dt,weight,csum,ssum;
  dt = getp("DTM",1,u);

  if(dt==0.0){
    csum=1.0/(2*N)*(1.0+u[3*(N-1)]);
    ssum=1.0/(2*N)*(0.0+u[3*(N-1)+N]);
    for (int k=0; k<N-1; k++){
      csum+=1.0/(2*N)*(u[k]+u[2*(N-1)+k]);
      ssum+=1.0/(2*N)*(u[(N-1)+k]+u[2*(N-1)+N+k]);
    }
    order=csum*csum+ssum*ssum;
    par[10]=0;
  }

  else{
    order=0;
    for (int i=0; i<NTST; i++){
      dt = getp("DTM",i+1,u);
      for (int j=0; j<NCOL+1; j++){
        weight = getp("WINT",j,u);
        csum=1.0/(2*N)*(1.0+ARRAY2D(u,3*(N-1),NCOL*i+j));
        ssum=1.0/(2*N)*(0.0+ARRAY2D(u,3*(N-1)+N,NCOL*i+j));
        for (int k=0; k<N-1; k++){
          csum+=1.0/(2*N)*(ARRAY2D(u,k,NCOL*i+j)+ARRAY2D(u,2*(N-1)+k,NCOL*i+j));
          ssum+=1.0/(2*N)*(ARRAY2D(u,(N-1)+k,NCOL*i+j)+ARRAY2D(u,2*(N-1)+N+k,NCOL*i+j));
        }
        order+=dt*weight*(csum*csum+ssum*ssum);
      }
    }
  }
  par[1]=pow(order,0.5);
  par[2]=getp("STA",0,u);
  par[3]=fabs(getp("STP",0,u));
  double mx=0,mr=0,mi=0;
  for(int i=1; i<ndim; i++){
    double re=getp("EIG",2*i+1,u);
    double im=getp("EIG",2*i+2,u);
    if(re*re+im*im>mx){
      mx=re*re+im*im;
      mr=re;
      mi=im;
    }
  }
  par[4]=mr;
  par[5]=mi;

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
