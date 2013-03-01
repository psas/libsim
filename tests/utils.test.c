#include <stdio.h>
#include <math.h>
#include "../libsim_types.h"
#include "test.h"
#include "../utils/coord.h"
#include "../physics/models/earth.h"
#include "utils.test.h"


char *ECEF2GEO_test(void)
{
	// Simple location
	vec ecef = {.v={ 0, RADIUS_EARTH, 0 }};
	vec geo = ECEF2GEO(ecef);
	double lon_err = fabs(geo.v.i - 0);
	double lat_err = fabs(geo.v.j - PI/2.0);
	double alt_err = fabs(geo.v.k - 0);

	printf("%f,%f,%f\n", geo.v.i, geo.v.j, geo.v.k);

	double tolerance = 1.0e-5;

	char * err = "\n  (-) Error: ECEF2GEO_test()\n        (+) Broken Transform\n";
	mu_assert(err, lon_err < tolerance);
	mu_assert(err, lat_err < tolerance);
	mu_assert(err, alt_err < 1);

	return 0; // tests passed
}
