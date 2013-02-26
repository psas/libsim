/**
 * Thrust 
 *
 * Author: Nathan Bergey
 * Version: Thursday, February 10 2011
 */
#include <stdio.h>
#include <stdlib.h> 
#include "../libsim_types.h"
#include "../math/vector.h"
#include "../utils/coord.h"
#include "models/earth.h"
#include "models/thrustcurve.h"
#include "thrust.h"


/**
 * Functions
 */
static double get_thrust_curve_segment(double t);

/**
* Thrust
*/
vec thrust(state s, double t, double *mdot)
{
  vec d;
  
  (*mdot) = get_thrust_curve_segment(t);
  
  double calc_thrust = rocket_thrust.Isp * g_0 * (*mdot);
  
  //printf("%f, %f, %f\n", t, thrust, *mdot);
  d.v.i = calc_thrust;
  d.v.j = 0;
  d.v.k = 0;
  
  return d;
}

void set_thrust_curve(thrust_curve calc_thrust)
{
  rocket_thrust = calc_thrust;
}

double get_thrust_curve_segment(double t)
{
  if (t < 0)
    return rocket_thrust.m_dot[0];
  else if (t > rocket_thrust.time[rocket_thrust.n])
    return 0;
  
  int i;
  for (i=0;i<=rocket_thrust.n;i++)
  {
    if (rocket_thrust.time[i] >= t)
    {
      //printf("%d: %f - %f  %f\n", i, rocket_thrust.time[i], t, rocket_thrust.m_dot[i]);
      break;
    }
  }

  return rocket_thrust.m_dot[i];  
}


void build_thrust_curve(double fuel, double isp, double avg_thrust, thrust_curve *t)
{
  double resolution = 0.1;

  double mdot = avg_thrust / (g_0*isp);

  int n = (fuel / mdot) / resolution;

  double burn_time = (fuel / mdot);
  
  (*t).time = (double *) malloc(sizeof(double)*n);
  (*t).m_dot = (double *) malloc(sizeof(double)*n);
  (*t).n = n;
  (*t).Isp = isp;
  int i;
  for (i=0;i<n;i++)
  {
    double time = i*resolution;
    (*t).time[i] = time;
    (*t).m_dot[i] = mdot;
  }
  
  (*t).time[n] = burn_time;
  (*t).m_dot[n] = mdot;
}
