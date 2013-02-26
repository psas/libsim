/**
 * Maximum number of steps one integration is allowed to take
 */
#define MAXSTEPS 10000

//void astraea(state init, rocket_stage stage, double x1, double x2, state_history *history, double *time_taken, state_change (*physics)(state,double), double (*set_boundary_condition)(state *history, double *x_history, int point_history, double x2));
state_history astraea(state init, rocket_stage stage, double x1, double x2, double *time_taken, integration_strategy s);
void interpolate_history(double *time, state *states, int steps, float dt, state *history);

/**
 * Integration Error Tolorance
 */
static const double eps = 1e-6;
