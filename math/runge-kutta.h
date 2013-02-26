#define SIGN(x,s) (((s)>0) ? (x) : (-(x)))      /* fortran definition */
#define FMAX(x,y)  ( (x)>(y) ? (x) : (y) )
#define FMIN(x,y)  ( (x)<(y) ? (x) : (y) )

#define TINY 1e-20
#define HSHRINK -0.25
#define HGROW -0.20
#define SAFETY 0.95
// TODO: Calc this
#define ERRCON 1.89e-4
//#define ERRCON 1e-10

void ode_int_fix_step(double *y, double *dydx, double *x, double h, int n, 
  void(*f)(double[],double[],double));
void rk4(double y[], double f1[], double x, int n, double h,
  void(*f)(double[],double[],double));
void rkqc(double *y, double *dydx, double *x, double htry, double eps, 
	double *yscal, double *hdid, double *hnext, int n,
	void(*f)(double[],double[],double));
