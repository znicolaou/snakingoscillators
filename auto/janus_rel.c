#include "auto_f2c.h"

int func (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
  /* System generated locals */
  // integer dfdu_dim1 = ndim;
  // integer dfdp_dim1 = ndim;

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

  // for(j=0; j<N-1; j++){
  //   printf("%f ",f[j]);
  // }
  // printf("\n");
  //
  // for(j=N-1; j<2*(N-1); j++){
  //   printf("%f ",f[j]);
  // }
  // printf("\n");
  //
  // for(j=2*(N-1); j<2*(N-1)+N; j++){
  //   printf("%f ",f[j]);
  // }
  // printf("\n");
  //
  // for(j=2*(N-1)+N; j<2*(N-1)+2*N; j++){
  //   printf("%f ",f[j]);
  // }
  // printf("\n");
  // exit(0);

  if (ijac == 0) {
    return 0;
  }
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
  // par[11] = 110;
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int pvls (integer ndim, const doublereal *u,
          doublereal *par)
{

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
