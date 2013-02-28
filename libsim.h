/**
 * Maximum number of steps one integration is allowed to take
 */
#define MAXSTEPS 10000

void Init_Model(void);
state_history Integrate_Rocket(rocket r, state initial_conditions);

/**
 * Integration Error Tolorance
 */
static const double eps = 1e-6;
