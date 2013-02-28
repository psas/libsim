/**
 * Physics
 *
 * Author: Nathan Bergey
 */

/**
 * Includes
 */
#include <stdio.h>
#include "../libsim_types.h"
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
