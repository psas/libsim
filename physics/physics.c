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
#include "../utils/structs.h"
#include "models/earth.h"
#include "gravity.h"
#include "aero.h"
#include "thrust.h"
#include "../libsim.h"
#include "physics.h"

/**
 * Functions
 */
state_change physics(state s, double t, integration_strategy strategy)
{
  state_change model;

  vec g = strategy.gravity_model(s);
 
  model.acc.v.i = (g.v.i) / s.m;
  model.acc.v.j = (g.v.j) / s.m;
  model.acc.v.k = (g.v.k) / s.m;
  model.m_dot = 0;
  return model;
}

/*
state_change physics_ThreeDOF_Sphere_G_only(state s, double t)
{
  state_change model;
  vec g = gravity(s);
  
  model.acc.i = g.i / s.m;
  model.acc.j = g.j / s.m;
  model.acc.k = g.k / s.m;
  model.m_dot = 0;
  
  return model;
}

state_change physics_OneDOF_Fixed_G_No_Atmosphere(state s, double t)
{
  state_change model;
  double g = g_0;
  
  model.acc.i = -g;
  model.acc.j = 0;
  model.acc.k = 0;
  model.m_dot = 0;
  
  return model;
}

state_change physics_OneDOF_Fixed_G_Exp_Atmosphere(state s, double t)
{
  state_change model;
  vec d = drag(s);
  double g = g_0;
  
  model.acc.i = (d.i / s.m) - g;
  model.acc.j = 0;
  model.acc.k = 0;
  model.m_dot = 0;
  
  return model;
}

state_change physics_OneDOF_Fixed_G_Thrust_No_Atmosphere(state s, double t)
{
  state_change model;
  double g = g_0;
  double mdot = 0;
  vec th = thrust(s, t, &mdot);

  model.acc.i = (th.i / s.m) - g;
  model.acc.j = 0;
  model.acc.k = 0;
  model.m_dot = -mdot;
  
  return model;
}

state_change physics_OneDOF_Thrust_Exp_Atmosphere(state s, double t)
{
  state_change model;
  vec g = gravity(s);
  vec d = drag(s);
  double mdot = 0;
  vec th = thrust(s, t, &mdot);

  model.acc.i = (g.i + d.i + th.i) / s.m;
  model.acc.j = 0;
  model.acc.k = 0;
  model.m_dot = -mdot;
  return model;
}

state_change physics_OneDOF_Thrust_No_Atmosphere(state s, double t)
{
  state_change model;
  vec g = gravity(s);
  double mdot = 0;
  vec th = thrust(s, t, &mdot);

  model.acc.i = (g.i + th.i) / s.m;
  model.acc.j = 0;
  model.acc.k = 0;
  model.m_dot = -mdot;
  
  return model;
}
*/
