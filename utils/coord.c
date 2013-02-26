#include <math.h>
#include "structs.h"
#include "../physics/models/earth.h"
#include "../math/vector.h"
#include "coord.h"

vec ECEF2GEO(vec v)
{
  double r   = 0;
  double lat = 0;
  double lon = 0;
  vec enu;
  r = norm(v);
  if (r > 0.0)
   lat = acos(v.v.k / r);
  if (v.v.i != 0.0 && v.v.j != 0)
   lon = atan2(v.v.j, v.v.i);

  enu.v.i = lon;
  enu.v.j = lat;
  enu.v.k = r - RADIUS_EARTH;
  
  return enu;
}

vec GEO2ECEF(vec v)
{
  double x = 0;
  double y = 0;
  double z = 0;
  vec ecef;
  double r = v.v.k + RADIUS_EARTH;
  
  double ci = cos(v.v.i);
  double sj = sin(v.v.j);
  double si = sin(v.v.i);
  double cj = cos(v.v.j);
  
  x = r * ci * sj;
  y = r * si * sj;
  z = r * cj;
  
  ecef.v.i = x;
  ecef.v.j = y;
  ecef.v.k = z;
  
  return ecef;
}

/**
 * ECEF vector to ENU vector
 * Depends on location on Earth, provided as lat/lon pair
 *
 * @param v     A vector in ECEF coordinates
 * @param lon   A vector in ECEF coordinates
 * @param lat   A vector in ECEF coordinates
 */
vec ECEF2ENU(vec v, double lon, double lat)
{
  mat3 T;
  double slon = sin(lon);
  double clon = cos(lon);
  double slat = sin(lat);
  double clat = cos(lat);
  
  vec enu;
  
  T.m.x1 =   -slon;    T.m.y1 =    clon;    T.m.z1 =   0;
  T.m.x2 = -slat*clon; T.m.y2 = -slat*slon; T.m.z2 = clat;
  T.m.x3 =  clat*clon; T.m.y3 =  clat*slon; T.m.z3 = slat;
  
  enu = matrix_mult(T, v);
  return enu;
}

double altitude(vec ecef_x)
{
  return norm(ecef_x) - RADIUS_EARTH;
}

double vertical_velocity(state r)
{
  vec velocity = r.v;
  vec unit_pos = unit_vec(r.x);

  double v_vel = dot_prod(unit_pos, velocity);
  
  return v_vel;
}

double vertical_acceleration_gee(vec a)
{
  double g0 = 9.81;
  double acc = a.v.i;
  
  return acc / g0;
}
