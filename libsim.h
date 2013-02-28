/**
 * Maximum number of steps one integration is allowed to take
 */
#define MAXSTEPS 10000

state Integrate(rocket r, state initial_conditions, state_history *history);

/**
 * Integration Error Tolorance
 */
static const double eps = 1e-6;
