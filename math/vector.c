/**
 * Vector Math
 *
 * Author: Nathan Bergey
 * Version: Thursday, February 10 2011
 */

/**
 * Includes
 */
#include <stdio.h>
#include <math.h>
#include "../utils/structs.h"
#include "vector.h"

/**
 * norm
 */
double norm(vec v)
{
  double n;
  n = sqrt( (v.v.i*v.v.i) + (v.v.j*v.v.j) + (v.v.k*v.v.k) );
  return n;
}

/**
 * unit_vec
 */
vec unit_vec(vec v)
{
  //int i;
  double magnitude;
  vec unitVector;
    
  magnitude = norm(v);
  
  if (magnitude == 0)
  {
    unitVector.v.i = 0.0;
    unitVector.v.j = 0.0;
    unitVector.v.k = 0.0;
    return unitVector;
  }
  unitVector.v.i = v.v.i / magnitude;
  unitVector.v.j = v.v.j / magnitude;
  unitVector.v.k = v.v.k / magnitude;
  
  return unitVector;
}

/**
 * Dot Product
 */
double dot_prod(vec a, vec b)
{
  double dot = 0;
  
  dot += a.v.i*b.v.i;
  dot += a.v.j*b.v.j;
  dot += a.v.k*b.v.k;
  
  return dot;
}

/**
 * Scale a vector; multiply a vector by a scalar
 */
vec vec_scale(vec v, double s)
{
  vec scaled;
  
  scaled.v.i = v.v.i * s;
  scaled.v.j = v.v.j * s;
  scaled.v.k = v.v.k * s;
  
  return scaled;
}
/**
 * Vector times Matrix
 */
vec matrix_mult(mat3 m, vec v)
{
  vec ans;
  int i, j;
  
  // Init
  for (i=0;i<3;i++) ans.component[i] = 0.0;
  
  for (i=0;i<3;i++)
  {
    for (j=0;j<3;j++)
    {
      ans.component[i] += m.component[i][j] * v.component[j];
    }
  }

  return ans;
}

/** 
 * Convert an axis angle to an equivalent rotation matrix
 *
 * Often it is useful to define a rotation around some arbitrary axis
 * Given a vector whos direction is the axis of rotation and whos magnitude is
 * the angle to rotate through this returns an equvilant rotation matrix
 *
 * see:
 * http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToMatrix/
 */
mat3 axis_angle_to_rotation_matrix(vec axis_angle)
{
	mat3 dst;
	double x = axis_angle.v.i, y = axis_angle.v.j, z = axis_angle.v.k;

	double angle = norm(axis_angle);
	if(fabs(angle) < 1e-30)
		return (mat3) {{ .x1 = 1, .y2 = 1, .z3 = 1 }};
	x /= angle;
	y /= angle;
	z /= angle;

	double c = cos(angle);
	double s = sin(angle);
	double t = 1.0 - c;

	dst.component[0][0] = c + x*x*t;
	dst.component[1][1] = c + y*y*t;
	dst.component[2][2] = c + z*z*t;

	double tmp1 = x*y*t;
	double tmp2 = z*s;
	dst.component[1][0] = tmp1 + tmp2;
	dst.component[0][1] = tmp1 - tmp2;
	tmp1 = x*z*t;
	tmp2 = y*s;
	dst.component[2][0] = tmp1 - tmp2;
	dst.component[0][2] = tmp1 + tmp2;
	tmp1 = y*z*t;
	tmp2 = x*s;
	dst.component[2][1] = tmp1 + tmp2;
	dst.component[1][2] = tmp1 - tmp2;

	return dst;
}

