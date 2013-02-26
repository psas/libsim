#include "../utils/structs.h"
#include "../utils/coord.h"
#include "../math/vector.h"
#include "test.h"
#include "misc.test.h"

/**
 * @test This tests a 90 degree CW turn around the Z axis
 */
char * coord_transform_test2()
{
  vec v = {0,1,0};
  vec axis = {0,0,-Pi/2.0}; // 90 deg CW about Z axis
  vec a;
  mat3 M;
  
  M = axis_angle_to_rotation_matrix(axis);
  a = matrix_mult(M, v);
  
  char * err = "\n  (-) Error: coord_transform_test2()\n        (+) Failed rotation\n";
  
  mu_assert(err, 
       (   (a.component[0] == 1 )
        && (a.component[1] < 1e-16 )
        && (a.component[1] < 1e-16 )));

  return 0;
}

/**
 * @test This tests the creation of a 90 degree rotation matrix, but does not
 * do the actual rotation.
 */
char * coord_transform_test1()
{
  vec v = {0,1,0};
  vec a;
  mat3 M;
  
  // 90deg CW
  M.x1 = 0; M.y1 = 1; M.z1 = 0;
  M.x2 = -1; M.y2 = 0; M.z2 = 0;
  M.x3 = 0; M.y3 = 0; M.z3 = 1;
  
  a = matrix_mult(M, v);

  char * err = "\n  (-) Error: coord_transform_test1()\n        (+) Failed rotation matrix creation\n";
  
   mu_assert(err, ( (a.component[0] == 1 ) 
                 && (a.component[1] == 0 ) 
                 && (a.component[2] == 0 )));

  return 0;
}
