/**
 * Main
 *
 * Author: Nathan Bergey
 */

#include <stdio.h>
#include <stdlib.h>
#include "libsim_types.h"
#include "libsim.h"

/**
 * main
 */
int main(int argc, char **argv)
{
	int i;

	/// Initilize libsim models
	Init_Model();

	/// Initilize a rocket
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
	state initial_conditions = { .x = {.v={-2414.59e3, -3771.092e3, 4528.117e3}},
                                 .v = {.v={0,0,0}},
                                 .a = {.v={0,0,0}},
                                 .m = 45
                               };
	state_history flight_history;

    /// Run a simulation
	flight_history = Integrate_Rocket(a_rocket, initial_conditions);

	for (i=0;i<flight_history.length;i++)
	{
		double x = flight_history.times[i];
		state  h = flight_history.states[i];
		printf("x: %f,    y:%f\n", x, h.x.v.i);
	}

	return 0; //exit
}
