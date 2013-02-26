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
 * @brief All of the Data Structures for the program
 *
 * @section DESCRIPTION
 *
 * Creates the data structures for Astraea
 */
 
/**
 * @brief Vector (3)
 *
 * A definition of a generic vector of length 3
 */
typedef union vec {
	struct {
		double i, j, k;
	} v;
	double component[3];
} vec;

/**
 * Matrix
 */
typedef union mat3
{
	struct {
		double x1, y1, z1;
		double x2, y2, z2;
		double x3, y3, z3;
	} m;
	double component[3][3];
} mat3;

/**
 * Rocket State
 */
typedef struct {vec x; vec v; vec a; double m;} state;

/**
 * Thrust Curve
 */
typedef struct {double *time; double *m_dot; int n; double Isp;} thrust_curve;

/**
 * Used to return the state from the physics model
 */
typedef struct {vec acc; double m_dot;} state_change;

/**
 * Used to return the an arrany of states and times from the integration
 */
typedef struct {double *times; state *states; int length;} state_history;

/**
 * Represents a single rocket stage
 */
typedef struct {thrust_curve thrust; double area; double Cd;} rocket_stage;

/**
 * A collection of rocket stages with common properties
 */
typedef struct {rocket_stage *stages; int number_of_stages;} rocket;


typedef void   (*integrator)(state *yp, double *xp, double x1, double x2, int *steps);
typedef vec    (*gravity)   (state s);
typedef vec    (*aero)      (state s);
typedef double (*boundary_condition)(state *history, double *x_history, int point_history, double x2);

/**
 *
 */
typedef struct {
  integrator integrator;
  int DOF;
  boundary_condition boundary;
  gravity gravity_model;
  aero    drag_model;
  } integration_strategy;


/**
 * PI
 */
#define Pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
