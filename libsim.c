/**
 * @file
 * @author  Nathan Bergey <nathan.bergey@gmail.com>
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
 * Calls to the rest of the program are handled through here.
 */

#include <stdio.h>
#include <stdlib.h>
#include "libsim_types.h"
#include "physics/physics.h"
#include "libsim.h"

/// Number of ODE's (3 position, 3 velocity, mass)
#define n 7

//static void state_to_rkvectors(double *y, double *dydx);
//static void rkvectors_to_state(double *y, double *dydx);
//void deriv(double *y ,double *dydx, double t);
//static void integrate(state *yp, double *xp, double x1, double x2, int *steps);

state Integrate(rocket r, state initial_conditions, state_history *history)
{
	/// Set the current state from the incoming initial state
	state current_state = initial_conditions;

	/// Return the final state after integrating
	return current_state;
}

/**
 * Integrate from end to end
 *
 * Does the integration
 *
 * @param yp An array of states to be apended to
 * @param xp An array of times to be appended to
 * @param x1 Begining time of the integration
 * @param x2 Ending time of the integration
 * @param steps Retunrs the number of steps the integrator took (does not
 *         include retried steps from rkck)
 */
/*
void integrate(state *yp, double *xp, double x1, double x2, int *steps)
{
	return;
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
*/
