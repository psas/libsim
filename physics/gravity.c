/**
 * Gravity
 *
 * Author: Nathan Bergey
 * Version: Thursday, February 10 2011
 */

/**
 * Includes
 */
#include <stdio.h>
#include "../utils/structs.h"
#include "../math/vector.h"
#include "models/earth.h"
#include "gravity.h"

/**
 * Gravity
 */
vec gravity_sphere(state s)
{
  vec g,e;
  double calc_gravity;
  double r = norm(s.x);
  
  // G Mm/r^2
  calc_gravity = G * ( (MASS_EARTH * s.m)/(r*r) );
  
  // Direction of gravity vector
  e = unit_vec(s.x);
  g.v.i = calc_gravity * e.v.i;
  g.v.j = calc_gravity * e.v.j;
  g.v.k = calc_gravity * e.v.k;
  
  return g;
}
