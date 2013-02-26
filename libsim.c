/**
 * @file
 * @author  Nathan Bergey <nathan.bergey@gmail.com>
 * @version Monday, May 02 2011
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
 * @brief Main program
 *
 * @section DESCRIPTION
 *
 * Calls to the rest of the program are handled through here. Specifically call
 * astraea() to start an integration. Call as many times as you need to optimize
 * a parameter of the state.
 */
 /**
 * @mainpage Astraea — A fast integrator/optimizer
 *
 * @author Nathan Bergey nathan.bergey@gmail.com
 *
 * @section intro Introduction
 * Astraea can integrate a system of Ordinary Differential Equations 
 * (ODEs) very fast using an adaptive stepsize Runge-Cutta numerical integrator.
 *
 * -# astraea()
 *
 * <hr>
 * @section requirements Requirements
 * @verbinclude requirements
 */

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>
#include "utils/structs.h"
#include "physics/physics.h"
#include "physics/thrust.h"
#include "physics/models/earth.h"
#include "math/runge-kutta.h"
#include "math/interpolation.h"
#include "utils/coord.h"
#include "libsim.h"

/// Current State of the rocket.
state current_state;

/// Number of ODE's (3 position, 3 velocity, mass)
#define n 7

static void state_to_rkvectors(double *y, double *dydx);
static void rkvectors_to_state(double *y, double *dydx);
void deriv(double *y ,double *dydx, double t);
static void integrate(state *yp, double *xp, double x1, double x2, int *steps);

integration_strategy strategy;

/**
 * Astraea main function
 * 
 * This is the entry point of the program, called by either main.c or can be
 * wrapped by other programs / languages
 *
 * @param init          The initial state of the rocket
 * @param stage         A rocket to simulate
 * @param x1            The begning time of the integration
 * @param x2            The ending time of the integration (if no other boundary conditions are met)
 * @param time_taken    Returns how much time the integration took
 */
state_history astraea(state init, rocket_stage stage, double x1, double x2, double *time_taken, integration_strategy s)
{
  state_history history;
  
  /// Set the current state from the incoming initial state
  current_state = init;

  strategy = s;
  //printf("set strategy\n");
  
  //printf("mass: %f\n", init.m);
  //printf("%f, %f, %f\n", init.x.i, init.x.j, init.x.k);

  ///TODO: fix with closure?
  /// If thrust curve exists, set in physics model
  if (stage.thrust.time != NULL)
    set_thrust_curve(stage.thrust);
    
  /// Measure the time at the beggining of the run
  //struct timeval t;
  //gettimeofday(&t, NULL);
  //long start = (t.tv_sec*1000000) + t.tv_usec;

  // Memory for integrator
  state yp[MAXSTEPS];
  double xp[MAXSTEPS];
  int steps_taken;
  
  /// Run the integration
  integrate(yp, xp, x1, x2, &steps_taken);
  
  /// Record state history
  // allocate enough memeory for the numer of steps
  history.states = malloc(sizeof(state) *(steps_taken));
  history.times  = malloc(sizeof(double)*(steps_taken));
  
  // For each step in steps_taken
  int i;
  for (i=0;i<steps_taken;i++)
  {
    history.states[i] = yp[i];
    history.times[i]  = xp[i];
  }
  history.length = steps_taken;
  
  /// Measure the time at the end of the run
  //gettimeofday(&t, NULL);
  //long end = (t.tv_sec*1000000) + t.tv_usec;
  
  /// Return the time it took
  //(*time_taken) = (end - start);
  
  return history;
}

/**
 * Integrate from end to end
 *
 * Does the integration
 *
 * \bug Setting the inition stepsize to the boundary length doesn't work.
 *
 * @param yp An array of states to be apended to
 * @param xp An array of times to be appended to
 * @param x1 Begining time of the integration
 * @param x2 Ending time of the integration
 * @param steps Retunrs the number of steps the integrator took (does not
 *         include retried steps from rkck)
 */
void integrate(state *yp, double *xp, double x1, double x2, int *steps)
{
  int i;
  int stepnum = 0;
	double yscale[n];
	double hdid;
	double hnext;
	double x = x1;
	//double h = x2-x1; \\TODO : Why doesn't this work?
	double h = 0.1;
	double time_to_stop = x2;
	
	// Init
	int nok = 0;
	int nbad = 0;
	double y[n];
	double dydx[n];
	state_to_rkvectors(y, dydx);
  
	// Start integration loop, through at most MAXSTEPS
	//for (stepnum = 0; stepnum <= MAXSTEPS; stepnum++)
	while (stepnum <= MAXSTEPS)
	{
	  // First RHS call
		deriv(y, dydx, x);
		
		// Store results for each step
	  xp[stepnum] = x;
	  yp[stepnum] = current_state;
	  
	  //printf("%f,%f,%f,%f\n", x, h, current_state.x.i, current_state.v.i);
	  /*
	  double exact_xi = -0.5*g_0*x*x + 5923*x + RADIUS_EARTH;
	  double exact_vi = 5923 - g_0*x;
	  printf("%f\n", g_0);
    printf("brons1,%f,%f,%f,%f\n", x, h, current_state.x.i, current_state.v.i);
    printf("brons2,%f,%f,%f,%f\n", x, h, exact_xi, exact_vi);
    printf("~~~~,%e,%e\n",current_state.x.i - exact_xi, current_state.v.i-exact_vi);
	  */
	
		// Y-scaling. Holds down fractional errors
		/// @todo Figure out proper Y-scaling
		//for (i=0;i<n;i++) yscale[i] = fabs(y[i]) + fabs(h*dydx[i]) + TINY;
		//for (i=0;i<n;i++) yscale[i] = fabs(h*dydx[i]) + TINY;
		for (i=0;i<n;i++) yscale[i] = dydx[i] + TINY;
		
	  //if (h > 0.5) h = 0.5;
			  
	  // Check for stepsize overshoot
		if ((x + h) > x2)
			h = x2 - x;
	  if ((x + h) > time_to_stop)
			h = time_to_stop - x;

		// One quality controled integrator step
		rkqc(y, dydx, &x, h, eps, yscale, &hdid, &hnext, n, deriv);
		
		// Was our timestep rational? if not record
		if (hdid==h)
    	nok++;
    else
    	nbad++;

    // Check boundary conditions, retry if overstepping:
    int q;
    state past[3];
    double x_past[3];

    // Create a recent history
    for (q=0;q<3;q++)
    {
      past[q] = yp[stepnum - q];
      x_past[q] = xp[stepnum - q];
    }
    
    // Evaluate boundary
    double new_x2 = strategy.boundary(past, x_past, 3, x2);
    
    // boundary violation
    if (new_x2 != x2)
    {
      // new attmpted stoping point
      x2 = new_x2;
      
      //printf("Step retry\n");
      //printf("RETRY_ %f  %f\n",x2, new_x2);
      // Retry last step
      stepnum = stepnum - 1;
      current_state = past[1];
      x = x_past[1];
      state_to_rkvectors(y, dydx);
    }
    else
    {
      // We might have set x2 too soon in the previous step, so to be safe we extend it.
      /// @todo: Fix this boundary conditions
      //printf("%f\n", h);
      x2 = x2;
    }
    
   	// Are we finished?
   	if ( (x2 - x) <= 0.0001 )
   	{
   	  deriv(y, dydx, x);
	    xp[stepnum] = x;
	    yp[stepnum] = current_state;
	    
	    // Done (boundary condition)
			break;
   	}
    if ( (time_to_stop - x) <= 0.0001 )
   	{
   	  deriv(y, dydx, x);
	    xp[stepnum] = x;
	    yp[stepnum] = current_state;

      // Done (time condion)
			break;
   	}

    // set timestep for next go around
   	h = hnext;
    stepnum++;
	}
	// How many steps did it take
	(*steps) = stepnum + 1;
}

void interpolate_history(double *time, state *states, int steps, float dt, state *history)
{
  int i, k;
  
  double **position;
  position = (double **) malloc(3 * sizeof(double *));
  for(i=0;i<3;i++)
  {
    position[i] = (double *) malloc(steps * sizeof(double));
    for (k=0;k<steps;k++)
    {
      position[i][k] = states[k].x.component[i];
    }
  }
  
  double **velocity;
  velocity = (double **) malloc(3 * sizeof(double *));
  for(i=0;i<3;i++)
  {
    velocity[i] = (double *) malloc(steps * sizeof(double));
    for (k=0;k<steps;k++)
    {
      velocity[i][k] = states[k].v.component[i];
    }
  }
  
  double **acceleration;
  acceleration = (double **) malloc(3 * sizeof(double *));
  for(i=0;i<3;i++)
  {
    acceleration[i] = (double *) malloc(steps * sizeof(double));
    for (k=0;k<steps;k++)
    {
      acceleration[i][k] = states[k].a.component[i];
    }
  }
  
  double **ypa;
  ypa = (double **) malloc(3 * sizeof(double *));
  for(i=0;i<3;i++)
  {
    ypa[i] = (double *) malloc(steps * sizeof(double));
    spline(time-1, position[i]-1, steps, velocity[i][0], velocity[i][steps], ypa[i]);
    //int l;
    //for (l = 0; l < steps; l ++)
    //  printf("%d    %f\n",l , ypa[i][l]);
  }
  
  double x = time[0];
  double xf = time[steps-1];
  int steps_to_take = ceil((xf-x)/dt);
  for (k=0;k<=steps_to_take;k++)
  {
    vec pos;
    for (i=0;i<3;i++)
    {
      double splinted;
      //splint(time,  position[i], ypa[i]-1, steps, x, &splinted);
      splint(time,  position[i], acceleration[i]-1, steps, x, &splinted);
      pos.component[i] = splinted;
    }
    history[k].x = pos;

    if ( (x + dt) > xf)
    {
      x = xf;
    }
    else
      x += dt;
  }
}

void state_to_rkvectors(double *y, double *dydx)
{
  y[0] = current_state.x.v.i;
  y[1] = current_state.v.v.i;
  
  y[2] = current_state.x.v.j;
  y[3] = current_state.v.v.j;
  
  y[4] = current_state.x.v.k;
  y[5] = current_state.v.v.k;
  
  dydx[1] = current_state.a.v.i;
  dydx[3] = current_state.a.v.j;
  dydx[5] = current_state.a.v.k;
  
  y[6] = current_state.m;
}

void rkvectors_to_state(double *y, double *dydx)
{
  current_state.x.v.i = y[0];
  current_state.v.v.i = y[1];
  
  current_state.x.v.j = y[2];
  current_state.v.v.j = y[3];
  
  current_state.x.v.k = y[4];
  current_state.v.v.k = y[5];
  
  current_state.a.v.i = dydx[1];
  current_state.a.v.j = dydx[3];
  current_state.a.v.k = dydx[5];
  
  current_state.m = y[6];
}

void deriv(double *y ,double *dydx, double t)
{
  // Velocity
  dydx[0] = y[1];
  dydx[2] = y[3];
  dydx[4] = y[5];
  
  // Acceleration
  rkvectors_to_state(y,dydx);
  state_change deriv_state = physics(current_state, t, strategy);
  current_state.a = deriv_state.acc;
  dydx[6] = deriv_state.m_dot;
  state_to_rkvectors(y, dydx);
}
