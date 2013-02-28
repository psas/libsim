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
 * @brief Runs Unit Tests
 *
 * @section DESCRIPTION
 *
 * Calls all the Unit Tests
 */
#include <stdio.h>
#include <stdlib.h>
#include "integrator.test.h"
#include "test.h"

int tests_run = 0;

static char * all_tests(void) {

	// Run Integrator Tests:
	mu_run_test(OneDOF_balistic_test1);

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

	return 0; //exit
}
