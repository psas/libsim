/**
 * @file
 * @author  Nathan Bergey <nathan.bergey@gmail.com>
 * @version Wednesday, April 27 2011
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @brief Runge-Kutta Integrator
 *
 * @section DESCRIPTION
 *
 * Based on Numerical Recipies in C
 * You are not expected to understand this.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "runge-kutta.h"

// Kash-Carp Tableau
static const double a[7] = {0,0,0.2,0.3,0.6,1.0,0.875};
														 // 0    			 1      2   			3 				4 5
static const double b[7][6] = {{0,    		 0,    0,  				0,				0,0} //0
															,{0,     		 0,    0,  				0,				0,0} //1
															,{0,  		 0.2,    0,  			  0,				0,0} //2
															,{0,		 0.075,0.225,  			  0, 				0,0} //3
															,{0,  		 0.3, -0.9,			  1.2,		    0,0} //4
															,{0,-11.0/54.0,  2.5,-70.0/27.0,35.0/27.0,0} //5
															,{0, 1631.0/ 55296.0
																 ,  175.0/   512.0
																 ,  575.0/ 13824.0
																 ,44275.0/110592.0
																 ,  253.0/  4096.0}};

static const double c[7]  = {0, 37.0/378.0
                              , 0
                              , 250.0/621.0
                              , 125.0/594.0
                              , 0
                              , 512.0/1771.0};
static const double dc[7] = {0, (37.0/378.0) - (2825.0/27648.0)
															, 0
															, 250.0/621.0 - 18575.0/48384.0
															, 125.0/594.0 - 13525.0/55296.0
															, -277.0/14336.0
															, 512.0/1771.0 - 0.25};

static void rkck(double *y, double *ak1, double x, double h, double *yout, double *yerr, int n, void(*f)(double[],double[],double));


void ode_int_fix_step(double *y, double *dydx, double *x, double h, int n, 
  void(*f)(double[],double[],double))
{
  f(y,dydx,(*x));
  rk4(y,dydx,(*x),n,h,f);
  (*x) += h;
}

/**
 * @brief A quality controled Runge-Kutta Integrator
 * 
 * RKQC stands for Runge-Kutta Quality Control.  There is a loop that runs 
 * inside here that trys an integration, and if the results are within the error
 * tolerance then it returns the new state vector and some results of the 
 * integration
 *
 * @param y The current RK state vector.
 * @param dydx An RK vector of state derivitaves
 * @param x The current time
 * @param htry The timestep to try
 * @param eps Error tolerance
 * @param yscal A scaling vector for the error tolorance
 * @param hdid A placeholder for the timestep that was actually used
 * @param hnext A suggested timestep for the next go around
 * @param n The number of elements in the RK vectors
 * @param f A function that will evaluate the derivative of the RK vectors
 */
void rkqc(double *y, double *dydx, double *x, double htry, double eps, 
	double *yscal, double *hdid, double *hnext, int n,
	void(*f)(double[],double[],double))
{
	int i;
	double errmax, htemp, h;
	//double xnew;
	double yerr[n];
	double ytemp[n];
	h = htry;

  /// Begin Loop
	for (;;)
	{
    /// Run one step
    rkck(y, dydx, *x, h, ytemp, yerr, n, f);

    /// Find the element with the highest error
    errmax = 0.0;
    for (i=0;i<n;i++) 
    {
      errmax = FMAX(errmax, fabs(yerr[i]/yscal[i]));
    }
    errmax /= eps;

    /// If within error tolorences then we're done
    if (errmax <= 1.0) break;

    /// If the error is too high scale h
    htemp = SAFETY * h * pow(errmax, HSHRINK);
    // But don't scale by more than a factor of 10 (beware of sign)
    h = (h >= 0.0) ? FMAX(htemp, 0.1*h) : FMIN(htemp, 0.1*h);

    /// Check for underflow
    //xnew = (*x) + h;
    //if (xnew == *x)
    if (h < TINY)
		{
      printf("Stepsize underflow\n");
      exit(1);
		}
	}/// Repeat
	
	/// Loop exited cleanly so we can increase timestep for next go-round
	if (errmax > ERRCON)
	{
		*hnext = SAFETY * h * pow(errmax, HGROW);
	}
	else
	{
		*hnext = 5.0 * h; 	// Limit to 5 times growth
	}
	  
  /// Update x (and hdid)
	*hdid = h;
	*x += h;
	
	/// Update values
	for (i=0;i<n;i++) y[i] = ytemp[i];
}

void rkck(double *y, double *ak1, double x, double h, double *yout, double *yerr, int n,
	void(*f)(double[],double[],double))
{
	int i;
	
	double ak2[n], ak3[n], ak4[n], ak5[n], ak6[n];
	double ytemp[n];
	
	/// Steps:
	
	// First Step
	for (i=0;i<n;i++)
	  ytemp[i] = y[i] + h*(b[2][1]*ak1[i]);
	  
	// Second Step
	(*f)(ytemp, ak2, x + h*a[2]);
	for (i=0;i<n;i++)
	  ytemp[i] = y[i] + h*(b[3][1]*ak1[i] + b[3][2]*ak2[i]);
	
	// Third Step
	(*f)(ytemp, ak3, x + h*a[3]);
	for (i=0;i<n;i++)
	  ytemp[i] = y[i] + h*(b[4][1]*ak1[i] + b[4][2]*ak2[i] + b[4][3]*ak3[i]);
	
	// Fourth Step
	(*f)(ytemp, ak4, x + h*a[4]);
	for (i=0;i<n;i++)
	  ytemp[i] = y[i] + h*(b[5][1]*ak1[i] + b[5][2]*ak2[i] + b[5][3]*ak3[i] + b[5][4]*ak4[i]);
	
	// Fifth Step
	(*f)(ytemp, ak5, x + h*a[5]);
	for (i=0;i<n;i++)
	  ytemp[i] = y[i] + h*(b[6][1]*ak1[i] + b[6][2]*ak2[i] + b[6][3]*ak3[i] + b[6][4]*ak4[i] + b[6][5]*ak5[i]);
	  
	// Sixth Step
  (*f)(ytemp, ak6, x + h*a[6]);
  
  /// Accumulate 4th order solution
  for (i=0;i<n;i++)
    yout[i] = y[i] + h*(c[1]*ak1[i] + c[3]*ak3[i] + c[4]*ak4[i] + c[6]*ak6[i]);
  
  /// Estimate the error as the difference between the 4th and 5th order solutions
  for (i=0;i<n;i++)
    yerr[i] = h*(dc[1]*ak1[i] + dc[3]*ak3[i] + dc[4]*ak4[i] + dc[5]*ak5[i] + dc[6]*ak6[i]);
}

void rk4(double y[], double f1[], double x, int n, double h,
  void(*f)(double[],double[],double))
{
  int i;
  double f2[n], f3[n], f4[n], tmp[n];
  double hh = h/2.0;
  double h6 = h/6.0;
  double xh = x + hh;
  
  // First Step
  for (i=0;i<n;i++) tmp[i] = y[i] + hh*f1[i];
  
  // Second Step
  (*f)(tmp, f2, xh);
  for (i=0;i<n;i++) tmp[i] = y[i] + hh*f2[i];
  
  // Third Step
  (*f)(tmp, f3, xh);
  for (i=0;i<n;i++) tmp[i] = y[i] + h*f3[i];
  
  // Fourth Step
  (*f)(tmp, f4, x+h);
  
  // Add Up
  for (i=0;i<n;i++)
    y[i] += h6*(f1[i] + 2.0*(f2[i] + f3[i]) + f4[i]);
}
