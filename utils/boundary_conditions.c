/**
 * @file
 * @author  Nathan Bergey <nathan.bergey@gmail.com>
 * @version Saturday, April 30 2011
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
 * @brief Defines types of boundary conditions
 *
 * @section DESCRIPTION
 *
 * The integration is asked to run between a beginning time and an ending time
 * however sometimes we would like to stop the integration at some in between 
 * point based on some condition, for example hitting the gound.
 *
 * astraea() takes a pointer to a function to find the boundary contitions for 
 * the integration. This file contains common usefull and reusable conditons in 
 * the required format.
 */
#include <stdio.h>
#include <math.h>
#include "structs.h"
#include "../math/interpolation.h"
#include "coord.h"
#include "boundary_conditions.h"

/**
 * @brief Ground impact condition
 * 
 * Use this condition for when you want to stop the integration when the vehicle 
 * hits the ground.
 *
 * @param history An array of the last few states
 * @param x_history A corisponding array of the last few times
 * @param point_history Number of elements in the arrays
 * @param x2 The current guess for boundary condition time
 *
 * @returns A guess for the boundary condition time 
 */
double boundary_condition_ground(state *history, double *x_history, int point_history, double x2)
{
    double current_altitude = altitude(history[0].x);

    // Condition
    if (current_altitude < 0)
    {
      double last_alt = altitude(history[1].x);
      x2 = linear_interpolate(last_alt, x_history[1], current_altitude, x_history[0], 0.0);
    }

    return x2;
}

/**
 * @brief Maximum height condition
 * 
 * Use this condition for when you have a simple balistic trajectory that reaches
 * a maximum height and then falls. It will allow the simulation to break at the 
 * maximum height.
 *
 * @param history An array of the last few states
 * @param x_history A corisponding array of the last few times
 * @param point_history Number of elements in the arrays
 * @param x2 The current guess for boundary condition time
 *
 * @returns A guess for the boundary condition time 
 */ 
double boundary_condition_max_alt(state *history, double *x_history, int point_history, double x2)
{
    double current_velocity = vertical_velocity(history[0]);

    // Condition
    if (current_velocity < 0)
    {
      double last_vel = vertical_velocity(history[1]);
      x2 = linear_interpolate(last_vel, x_history[1], current_velocity, x_history[0], 0.0);
    }

    return x2;
}
