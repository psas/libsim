/**
 * Aero
 *
 * Author: Nathan Bergey
 * Version: Wednesday, April 13 2011
 */
#include <stdio.h>
#include <math.h>
#include "../libsim_types.h"
#include "../math/vector.h"
#include "../utils/coord.h"
#include "aero.h"

/**
* Drag
*/
vec drag(state s)
{
  vec d;
  double v = norm(s.v);
  vec v_hat = unit_vec(s.v);
  double h = altitude(s.x);
  double rho_0 = 1.2041;
  double k = 0.000115; // A guess
  double Cd = 1;
  double A = 0.015;
  double calc_drag = 0;
  double rho = rho_0;
  
  rho = rho_0*exp(-k*h);
  
  calc_drag = -0.5*rho*v*v*A*Cd;
  
  d = vec_scale(v_hat, calc_drag);
  
  return d;
}
