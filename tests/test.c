/**
 * @file
 * @author  Nathan Bergey <nathan.bergey@gmail.com>
 * @version Sunday, May 01 2011
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
 * @brief Runs Unit Tests
 *
 * @section DESCRIPTION
 *
 * Calls all the Unit Tests
 */
#include <stdio.h>
#include "integrator.test.h"
#include "misc.test.h"
#include "test.h"

int tests_run = 0;

static char * all_tests() {
  /// Coordinate Transform Tests
  mu_run_test(coord_transform_test1);
  mu_run_test(coord_transform_test2);

  /// Integrator Tests
  mu_run_test(balistic_test1);
  mu_run_test(thrust_test1);
  mu_run_test(thrust_test2);
  return 0;
}

/**
 * @brief Run
 */
int main(int argc, char **argv)
{
  char *result = all_tests();
  if (result != 0) {
    printf("%s\n", result);
    exit(1);
  }
  else {
    printf("\nALL %d TESTS PASSED\n\n", tests_run);
  }

  return 0;
}
