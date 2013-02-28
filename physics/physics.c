/**
 * Physics
 *
 * Author: Nathan Bergey
 * Version: Thursday, February 10 2011
 */

/**
 * Includes
 */
#include <stdio.h>
#include "../libsim_types.h"
//#include "models/earth.h"
#include "gravity.h"
#include "aero.h"
#include "thrust.h"
//#include "../libsim.h"
#include "physics.h"

/**
 * Functions
 */
state_change physics(state s, double t, physics_model_strategy strategy)
{
    // Return value
    state_change model;

    // Calc gravity
    vec g = strategy.gravity_model(s);

    model.acc.v.i = (g.v.i) / s.m;
    model.acc.v.j = (g.v.j) / s.m;
    model.acc.v.k = (g.v.k) / s.m;
    model.m_dot = 0;

    return model;
}
