#include "auto_f2c.h"

void quicksort(double *vec,int first,int last);

int func (integer ndim, const doublereal *w, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
  /* System generated locals */
  integer dfdu_dim1 = ndim;
  integer dfdp_dim1 = ndim;

  double amp = par[0];
  double omega = par[1];
  double delta= par[2];
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
    f[0*N+j]=-p[j]*y[j]/(1+delta)/omega+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
    f[1*N+j]= p[j]*x[j]/(1+delta)/omega+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
    f[2*N+j]=(-eta*p[j]-(1+omega*omega*amp*X+4*delta)*y[j]+(1-delta)*(x[j]*v[j]-u[j]*y[j]+x[j]*v[(j-1+N)%N]-y[j]*u[(j-1+N)%N]))/omega;
    f[3*N+j]=-q[j]*v[j]/(1-delta)/omega+gamma*(1-u[j]*u[j]-v[j]*v[j])*u[j];
    f[4*N+j]= q[j]*u[j]/(1-delta)/omega+gamma*(1-u[j]*u[j]-v[j]*v[j])*v[j];
    f[5*N+j]=(-eta*q[j]-(1+omega*omega*amp*X-4*delta)*v[j]+(1+delta)*(u[j]*y[j]-x[j]*v[j]+u[j]*y[(j+1)%N]-v[j]*x[(j+1)%N]))/omega;
  }
  f[6*N+0]=-Y+gamma*(1-X*X-Y*Y)*X;
  f[6*N+1]= X+gamma*(1-X*X-Y*Y)*Y;

  if (ijac == 0) {
    return 0;
  }

  for(j=0; j<N; j++){
    // f[0*N+j]=-p[j]*y[j]/(1+delta)/omega+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
    ARRAY2D(dfdu,0*N+j,0*N+j)=gamma*(1-3*x[j]*x[j]-y[j]*y[j]);
    ARRAY2D(dfdu,0*N+j,1*N+j)=-p[j]/(1+delta)/omega+gamma*(-2*y[j])*x[j];
    ARRAY2D(dfdu,0*N+j,2*N+j)=-y[j]/(1+delta)/omega;

    // f[1*N+j]= p[j]*x[j]/(1+delta)/omega+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
    ARRAY2D(dfdu,1*N+j,0*N+j)=p[j]/(1+delta)/omega+gamma*(-2*x[j])*y[j];
    ARRAY2D(dfdu,1*N+j,1*N+j)=gamma*(1-x[j]*x[j]-3*y[j]*y[j]);
    ARRAY2D(dfdu,1*N+j,2*N+j)=x[j]/(1+delta)/omega;

    // f[2*N+j]=(-eta*p[j]-(1+omega*omega*amp*X+4*delta)*y[j]+(1-delta)*(x[j]*v[j]-u[j]*y[j]+x[j]*v[(j-1+N)%N]-y[j]*u[(j-1+N)%N]))/omega;
    ARRAY2D(dfdu,2*N+j,0*N+j)=((1-delta)*(v[j]+v[(j-1+N)%N]))/omega;
    ARRAY2D(dfdu,2*N+j,1*N+j)=(-(1+omega*omega*amp*X+4*delta)+(1-delta)*(-u[j]-u[(j-1+N)%N]))/omega;
    ARRAY2D(dfdu,2*N+j,2*N+j)=(-eta)/omega;
    ARRAY2D(dfdu,2*N+j,3*N+j)=((1-delta)*(-y[j]))/omega;
    ARRAY2D(dfdu,2*N+j,4*N+j)=((1-delta)*(x[j]))/omega;
    ARRAY2D(dfdu,2*N+j,3*N+(j-1+N)%N)=((1-delta)*(-y[j]))/omega;
    ARRAY2D(dfdu,2*N+j,4*N+(j-1+N)%N)=((1-delta)*(x[j]))/omega;
    ARRAY2D(dfdu,2*N+j,6*N+0)=(-(omega*omega*amp)*y[j])/omega;

    // f[3*N+j]=-q[j]*v[j]/(1-delta)/omega+gamma*(1-u[j]*u[j]-v[j]*v[j])*u[j];
    ARRAY2D(dfdu,3*N+j,3*N+j)=gamma*(1-3*u[j]*u[j]-v[j]*v[j]);
    ARRAY2D(dfdu,3*N+j,4*N+j)=-q[j]/(1-delta)/omega+gamma*(-2*v[j])*u[j];
    ARRAY2D(dfdu,3*N+j,5*N+j)=-v[j]/(1-delta)/omega;

    // f[4*N+j]= q[j]*u[j]/(1-delta)/omega+gamma*(1-u[j]*u[j]-v[j]*v[j])*v[j];
    ARRAY2D(dfdu,4*N+j,3*N+j)=q[j]/(1-delta)/omega+gamma*(-2*u[j])*v[j];
    ARRAY2D(dfdu,4*N+j,4*N+j)=gamma*(1-u[j]*u[j]-3*v[j]*v[j]);
    ARRAY2D(dfdu,4*N+j,5*N+j)=u[j]/(1-delta)/omega;

    // f[5*N+j]=(-eta*q[j]-(1+omega*omega*amp*X-4*delta)*v[j]+(1+delta)*(u[j]*y[j]-x[j]*v[j]+u[j]*y[(j+1)%N]-v[j]*x[(j+1)%N]))/omega;
    ARRAY2D(dfdu,5*N+j,0*N+j)=((1+delta)*(-v[j]))/omega;
    ARRAY2D(dfdu,5*N+j,1*N+j)=((1+delta)*(u[j]))/omega;
    ARRAY2D(dfdu,5*N+j,0*N+(j+1)%N)=((1+delta)*(-v[j]))/omega;
    ARRAY2D(dfdu,5*N+j,1*N+(j+1)%N)=((1+delta)*(u[j]))/omega;
    ARRAY2D(dfdu,5*N+j,3*N+j)=((1+delta)*(y[j]+y[(j+1)%N]))/omega;
    ARRAY2D(dfdu,5*N+j,4*N+j)=(-(1+omega*omega*amp*X-4*delta)+(1+delta)*(-x[j]-x[(j+1)%N]))/omega;
    ARRAY2D(dfdu,5*N+j,5*N+j)=(-eta)/omega;
    ARRAY2D(dfdu,5*N+j,6*N+0)=(-(omega*omega*amp)*v[j])/omega;

  }

  // f[6*N+0]=-Y+gamma*(1-X*X-Y*Y)*X;
  ARRAY2D(dfdu,6*N+0,6*N+0)=gamma*(1-3*X*X-Y*Y);
  ARRAY2D(dfdu,6*N+0,6*N+1)=-1+gamma*(-2*Y)*X;

  // f[6*N+1]= X+gamma*(1-X*X-Y*Y)*Y;
  ARRAY2D(dfdu,6*N+1,6*N+0)=1+gamma*(-2*X)*Y;
  ARRAY2D(dfdu,6*N+1,6*N+1)=gamma*(1-X*X-3*Y*Y);


  if (ijac == 1) {
    return 0;
  }

  for(j=0; j<N; j++){
    // f[0*N+j]=-p[j]*y[j]/(1+delta)/omega+gamma*(1-x[j]*x[j]-y[j]*y[j])*x[j];
    ARRAY2D(dfdp,0*N+j,1)=p[j]*y[j]/(1+delta)/omega/omega;
    ARRAY2D(dfdp,0*N+j,2)=p[j]*y[j]/(1+delta)/omega/(1+delta);

    // f[1*N+j]= p[j]*x[j]/(1+delta)/omega+gamma*(1-x[j]*x[j]-y[j]*y[j])*y[j];
    ARRAY2D(dfdp,1*N+j,1)=-p[j]*x[j]/(1+delta)/omega/omega;
    ARRAY2D(dfdp,1*N+j,2)=-p[j]*x[j]/(1+delta)/omega/(1+delta);

    // f[2*N+j]=(-eta*p[j]-(1+omega*omega*amp*X+4*delta)*y[j]+(1-delta)*(x[j]*v[j]-u[j]*y[j]+x[j]*v[(j-1+N)%N]-y[j]*u[(j-1+N)%N]))/omega;
    ARRAY2D(dfdp,2*N+j,0)=(-(omega*omega*X)*y[j])/omega;
    ARRAY2D(dfdp,2*N+j,1)=(eta*p[j]-(-1+omega*amp*X-4*delta)*y[j]-(1-delta)*(x[j]*v[j]-u[j]*y[j]+x[j]*v[(j-1+N)%N]-y[j]*u[(j-1+N)%N]))/omega/omega;
    ARRAY2D(dfdp,2*N+j,2)=(-(4)*y[j]+(-1)*(x[j]*v[j]-u[j]*y[j]+x[j]*v[(j-1+N)%N]-y[j]*u[(j-1+N)%N]))/omega;

    // f[3*N+j]=-q[j]*v[j]/(1-delta)/omega+gamma*(1-u[j]*u[j]-v[j]*v[j])*u[j];
    ARRAY2D(dfdp,3*N+j,1)=q[j]*v[j]/(1-delta)/omega/omega;
    ARRAY2D(dfdp,3*N+j,2)=-q[j]*v[j]/(1-delta)/omega/(1-delta);

    // f[4*N+j]= q[j]*u[j]/(1-delta)/omega+gamma*(1-u[j]*u[j]-v[j]*v[j])*v[j];
    ARRAY2D(dfdp,4*N+j,1)=-q[j]*u[j]/(1-delta)/omega/omega;
    ARRAY2D(dfdp,4*N+j,2)=q[j]*u[j]/(1-delta)/omega/(1-delta);

    // f[5*N+j]=(-eta*q[j]-(1+omega*omega*amp*X-4*delta)*v[j]+(1+delta)*(u[j]*y[j]-x[j]*v[j]+u[j]*y[(j+1)%N]-v[j]*x[(j+1)%N]))/omega;
    ARRAY2D(dfdp,5*N+j,0)=(-(omega*omega*X)*v[j])/omega;
    ARRAY2D(dfdp,5*N+j,1)=(eta*q[j]-(-1+omega*amp*X+4*delta)*v[j]-(1+delta)*(u[j]*y[j]-x[j]*v[j]+u[j]*y[(j+1)%N]-v[j]*x[(j+1)%N]))/omega/omega;
    ARRAY2D(dfdp,5*N+j,2)=(-(-4)*v[j]+(1)*(u[j]*y[j]-x[j]*v[j]+u[j]*y[(j+1)%N]-v[j]*x[(j+1)%N]))/omega;
  }

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par)
{
  par[0] = 0.06;
  par[1] = 3.4;
  par[2]=0.25;
  int j, N = (ndim-2)/6;
  for (j=0; j<N; j++){
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
  doublereal dt,weight,norm1=0,norm2=0;
  double *vec=malloc(ndim*sizeof(double));
  doublereal det1=1.0,det2=1.0;


  dt = getp("DTM",1,u);

  if(dt==0.0){
    par[10]=0;
  }

  else{

    for (int i=0; i<NTST; i++){
      dt = getp("DTM",i+1,u);
      for (int j=0; j<NCOL+1; j++){
        weight = getp("WINT",j,u);
        for(int k=0; k<N; k++){
          norm1+=dt*weight*(ARRAY2D(u,N+k,NCOL*i+j)*ARRAY2D(u,N+k,NCOL*i+j)+ARRAY2D(u,4*N+k,NCOL*i+j)*ARRAY2D(u,4*N+k,NCOL*i+j));
          norm2+=dt*weight*(ARRAY2D(u,2*N+k,NCOL*i+j)*ARRAY2D(u,2*N+k,NCOL*i+j)+ARRAY2D(u,5*N+k,NCOL*i+j)*ARRAY2D(u,5*N+k,NCOL*i+j));
        }
      }
    }
    for(int i=0; i<ndim; i++){
      // vec[i]=log(sqrt(getp("EIG",2*i+1,u)*getp("EIG",2*i+1,u)+getp("EIG",2*i+2,u)*getp("EIG",2*i+2,u)));
      vec[i]=(getp("EIG",2*i+1,u)*getp("EIG",2*i+1,u)+getp("EIG",2*i+2,u)*getp("EIG",2*i+2,u)-1)/(getp("EIG",2*i+1,u)*getp("EIG",2*i+1,u)+getp("EIG",2*i+2,u)*getp("EIG",2*i+2,u)+1);
    }
    //Do an argsort so we can get real and imaginary parts of each sorted exponent
    quicksort(vec,1,ndim-1);

    det1=vec[1];
    det2=vec[2];
    for(int i=1; i<ndim; i++){
      if(vec[i]<0&&i!=1){
        det1=det1;
      }
      if(vec[i]<0&&i>2){
        det2=det2; //Fine unless vec[i] changes sign for i!=2...
      }
    }
  }

  // par[3]=det*(1+vec[2]);
  par[3]=getp("STA",0,u);
  par[4]=sqrt(norm1);
  par[5]=sqrt(norm2);
  par[6]=det1;
  par[7]=det2;

  free(vec);
  return 0;
}

void argsort(doublereal *vec, int *order, int ndim){
  double min=vec[0];
  //find the smallest element
  for (int i=1; i<ndim; i++){
    if(vec[i]<min){
      min=vec[i];
      order[0]=i;
    }
  }

  //find the next smallest value
  for (int i=1; i<ndim; i++){
    //find the smallest element of vector larger than the current min
      int ind=0;
      for(int j=0; j<ndim; j++){
        if(vec[j]<min && vec[j]>vec[order[i-1]]){
          min=vec[j];
          order[i]=j;
        }
      }
  }
}

void quicksort(doublereal *vec, int first, int last){
  int i, j, pivot;
  doublereal temp;

  if(first<last){
    pivot=first;
    i=first;
    j=last;

    while(i<j){
      // while(vec[i]<=vec[pivot]&&i<last){
      while(fabs(vec[i])<=fabs(vec[pivot])&&i<last){
        i++;
      }
      // while(vec[j]>vec[pivot]){
      while(fabs(vec[j])>fabs(vec[pivot])){
        j--;
      }
      if(i<j){
        temp=vec[i];
        vec[i]=vec[j];
        vec[j]=temp;
      }
    }

    temp=vec[pivot];
    vec[pivot]=vec[j];
    vec[j]=temp;
    quicksort(vec,first,j-1);
    quicksort(vec,j+1,last);
  }
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
