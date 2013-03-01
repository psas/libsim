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
#include <stdbool.h>
#include "libsim_types.h"
#include "physics/physics.h"
#include "physics/gravity.h"
#include "math/runge-kutta.h"
#include "libsim.h"

/// Number of ODE's (3 position, 3 velocity, mass)
#define n 7

// Local functions
void deriv(double *y ,double *dydx, double t);
static state rk2state(double *y, double *dydx);
static void integrate(state y0, state *yp, double *xp, double x1, double x2, int *steps);

// Globals
physics_model_strategy physics_model;

void Init_Model(void)
{
	physics_model.gravity_model = gravity_sphere;
}

state_history Integrate_Rocket(rocket r, state initial_conditions)
{
	int i;
	// History
	state *yp;
	double *xp;
	yp = malloc(sizeof(state) * MAXSTEPS);
	xp = malloc(sizeof(double) * MAXSTEPS);
	state_history history;
	state *h;
	double *times;

	int steps_taken = 0;

	integrate(initial_conditions, yp, xp, 0, 1, &steps_taken);

	h = malloc(sizeof(state) * steps_taken);
	times = malloc(sizeof(double) * steps_taken);
	for (i=0;i<steps_taken;i++)
	{
		h[i] = yp[i];
		times[i] = xp[i];
	}

	history.times = times;
	history.states = h;
	history.length = steps_taken;

	return history;
}

static void integrate(state y0, state *yp, double *xp, double x1, double x2, int *steps)
{
	int i;
	int stepnum = 0;   // Track number of steps the integrator has run

	double x = x1;     // Current time, begin at x1
	state s;           // Current state

	// Integrator memory, each position in an array is a DOF of the system
	double y[n];       // array of integrator outputs, y = integral(y' dx)
	double dydx[n];    // array of RHS, dy/dx
	double yscale[n];  // array of yscale factors (integraion error tolorence)
	double h;          // timestep
	double hdid;       // stores actual timestep taken by RK45
	double hnext;      // guess for next timestep

	// stop the integrator
	double time_to_stop = x2;

	// First guess for timestep
	h = x2 - x1;

	// Inital conditions
	y[0] = y0.x.v.i;
	y[1] = y0.v.v.i;
	y[2] = y0.x.v.j;
	y[3] = y0.v.v.j;
	y[4] = y0.x.v.k;
	y[5] = y0.v.v.k;
	dydx[1] = y0.a.v.i;
	dydx[3] = y0.a.v.j;
	dydx[5] = y0.a.v.k;
	y[6] = y0.m;

	while (stepnum <= MAXSTEPS)
	{
		// First RHS call
		deriv(y, dydx, x);
		s = rk2state(y, dydx);

		// Store current state
		xp[stepnum] = x;
		yp[stepnum] = s;

		// Y-scaling. Holds down fractional errors
		for (i=0;i<n;i++) {
			yscale[i] = 0.000000001;//dydx[i]+TINY;
		}

		// Check for stepsize overshoot
		if ((x + h) > time_to_stop)
			h = time_to_stop - x;

		// One quality controled integrator step
		rkqc(y, dydx, &x, h, eps, yscale, &hdid, &hnext, n, deriv);
		s = rk2state(y, dydx);

		// Are we finished?
		// hit ground
		if (underground(s))
		{
			deriv(y, dydx, x);
			s = rk2state(y, dydx);
			xp[stepnum] = x;
			yp[stepnum] = s;
			stepnum++;
			break;
		}

		// Passed requested integration time
		if ( (time_to_stop - x) <= 0.0001 )
		{
			deriv(y, dydx, x);
			s = rk2state(y, dydx);
			xp[stepnum] = x;
			yp[stepnum] = s;
			stepnum++;
			break;
		}

		// set timestep for next go around
		h = hnext;

		// Incement step counter
		stepnum++;
	}

	(*steps) = stepnum;
	return;
}

static state rk2state(double *y, double *dydx)
{
	state s;
	s.x.v.i = y[0];
	s.v.v.i = y[1];
	s.x.v.j = y[2];
	s.v.v.j = y[3];
	s.x.v.k = y[4];
	s.v.v.k = y[5];
	s.m     = y[6];

	s.a.v.i = dydx[1];
	s.a.v.j = dydx[3];
	s.a.v.k = dydx[5];

	return s;
}

void deriv(double *y ,double *dydx, double t)
{
	state current_state;

	// To solve for acceleration build a state struct to send to physics engine
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

	// Do Physics to current state:
	state_change deriv_state = physics(current_state, t, physics_model);

	// Build RK vectors from state change
	// Velocity is single integration of acceleration
	dydx[0] = y[1];
	dydx[2] = y[3];
	dydx[4] = y[5];
	// Acceleration is from physics model
	dydx[1] = deriv_state.acc.v.i;
	dydx[3] = deriv_state.acc.v.j;
	dydx[5] = deriv_state.acc.v.k;
	dydx[6] = deriv_state.m_dot;
}
