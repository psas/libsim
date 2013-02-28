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
	rocket initial_rocket = { .thrust = motor,
                              .area = 0.4,
                              .Cd = 0.8
                            };
	state initial_conditions = { .x = {.v={0,0,0}},
                                 .v = {.v={0,0,0}},
                                 .a = {.v={0,0,0}},
                                 .m = 45
                               };
	state final_state;
	state_history *flight_history;

    /// Run a simulation
	final_state = Integrate_Rocket(initial_rocket, initial_conditions, flight_history);

	return 0; //exit
}
