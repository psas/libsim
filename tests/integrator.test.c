#include <stdio.h>
#include <stdlib.h>
#include "../libsim_types.h"
#include "../libsim.h"
#include "../physics/models/earth.h"
#include "../utils/coord.h"
#include "test.h"
#include "integrator.test.h"

/**
 * @test This sets up a 1D balistics model and compairs the result from astraea()
 * to the true, analytical solution at all returned points.
 */
char *OneDOF_balistic_test1(void)
{
	int i;
	double intial_velocity = 343.94;
	double inital_height = 0.0;

	// Internal Exact value calculators
	double exact(double t)
	{
		double height = -0.5*g_0*t*t + intial_velocity*t + inital_height;
		return height;
	}
	double exact_v(double t)
	{
		double velocity = -g_0*t + intial_velocity;
		return velocity;
	}

	// Initilize libsim models
	Init_Model();

	// Initilize a rocket
	double t[2] = {0,1};
	double m[2] = {1,1};
	thrust_curve motor = { .time = t,
                           .m_dot = m,
                           .length = 2,
                           .Isp = 254
                         };
	rocket a_rocket = { .thrust = motor,
                              .area = 0.4,
                              .Cd = 0.8
                            };
	vec position = {.v={-2.14031, 0.79412, inital_height}};
	position = GEO2ECEF(position);
	vec velocity = {.v={0,0,intial_velocity}};
	velocity = ENU2ECEF(velocity, -2.14031, 0.79412);
	state initial_conditions = { .x = position,
                                 .v = {.v={0,0,0}},
                                 .a = {.v={0,0,0}},
                                 .m = 45
                               };

    /// Run a simulation
	state_history history = Integrate_Rocket(a_rocket, initial_conditions);

	for (i=0;i<history.length;i++)
	{
		double x_error = altitude(history.states[i].x) - exact(history.times[i]);
		//double v_error = history.states[k].v.v.i - exact_v(history.times[k]);

		char * err = "\n  (-) Error: OneDOF_balistic_test1()\n        (+) Integration error larger than expected\n";

		mu_assert(err, x_error < 0.001); //To the mm
		//mu_assert(err, v_error < 0.001); //To the mm/s
	}

	return 0; // tests passed
}
