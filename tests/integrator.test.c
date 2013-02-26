#include <stdio.h>
#include <math.h>
#include "../utils/structs.h"
#include "../utils/coord.h"
#include "../utils/boundary_conditions.h"
#include "../math/runge-kutta.h"
#include "../math/interpolation.h"
#include "../math/vector.h"
#include "../physics/models/earth.h"
#include "../physics/physics.h"
#include "../physics/thrust.h"
#include "../astraea.h"
#include "test.h"
#include "integrator.test.h"

state init;
state_history history;
rocket the_rocket;

void init_rocket_setup();

/**
 * @test This sets up a 1D balistics model and compairs the result from astraea()
 * to the true, analytical solution at all returned points.
 */
char * balistic_test1()
{
  double intial_velocity = 343.94;
  double inital_height = 0.0;
  
  init_rocket_setup();
  
  // Exact
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
 
  double sim_time;
  
  init.x.i += inital_height;
  init.v.i = intial_velocity;
  init.m = 25.5;
  
  rocket_stage sustainer;
  
  astraea(  init, 
            sustainer,
            0, 
            1000,
            &history,
            &sim_time,
            physics_OneDOF_Fixed_G_Exp_Atmosphere,
            boundary_condition_max_alt);

  int k;
  for (k=0;k<=history.length;k++)
  {
    double x_error = altitude(history.states[k].x) - exact(history.times[k]);
    double v_error = history.states[k].v.i - exact_v(history.times[k]);

    char * err = "\n  (-) Error: balistic_test1()\n        (+) Integration error larger than expected\n";
    
    mu_assert(err, x_error < 0.001); //To the mm
    mu_assert(err, v_error < 0.001); //To the mm/s
  }
  
  return 0;
}

/**
 * @test This sets up a 1D thrust model and compairs the result 
 * from astraea() to the true, analytical solution at all returned points.
 */
char * thrust_test1()
{ 
  init_rocket_setup();
   
  double intial_velocity = 0.0;
  double inital_height = 0.0;
  double Isp = 260;
  double m_0 = 25.4;
  double avg_thrust = 2500;
  double fuel = 4.5;
  
  // Exact
  double exact_x(double t)
  {
    double mdot = avg_thrust/(g_0*Isp);
    double mass = m_0 - mdot * t;
    
    double alt = g_0*Isp*(m_0/mdot)*(1-(mass/m_0)*(log(m_0/mass)+1)) - 0.5*g_0*t*t;
    
    return alt;
  }
  
  double exact_v(double t)
  {
    double mass = m_0 - avg_thrust/(g_0*Isp) * t;
    double velocity = g_0*Isp*log(m_0 / mass) - g_0*t + intial_velocity;
    return velocity;
  }
  
  double exact_m(double t)
  {
    double mass = m_0 - avg_thrust/(g_0*Isp) * t;
    return mass;
  }
  
  double exact_a(double t)
  {
    double mdot = avg_thrust/(g_0*Isp);
    double mass = m_0 - mdot * t;
    double acc = ((g_0*Isp)/mass) * mdot - g_0;
    return acc;
  }
 
  double exact_burn_time = fuel / (avg_thrust/(g_0*Isp));
  
  double sim_time;
  
  init.x.i += inital_height;
  init.v.i = intial_velocity;
  init.m = m_0;
  
  rocket_stage sustainer;
  sustainer.area = 0.1;
  sustainer.Cd = 0.8;
  
  thrust_curve th;
  build_thrust_curve(fuel, Isp, avg_thrust, &th);
  
  double burn_time = th.time[th.n];
 
  //printf("%f == %f\n\n",burn_time, exact_burn_time);
  char * err = "\n  (-) Error: thrust_test1()\n        (+) Incorrect Burn Time\n";
  double burn_time_error = fabs(burn_time - exact_burn_time);
  mu_assert(err, burn_time_error < 0.00001);
  
  sustainer.thrust = th;
  
  astraea(  init, 
            sustainer,
            0,
            burn_time,
            &history,
            &sim_time,
            physics_OneDOF_Fixed_G_Thrust_No_Atmosphere,
            boundary_condition_max_alt);

  int k;
  for (k=0;k<=history.length-1;k++)
  {
    double t =  history.times[k];
    
    double x_error = fabs(altitude(history.states[k].x) - exact_x(t));
    double v_error = fabs(history.states[k].v.i - exact_v(t));
    double m_error = fabs(history.states[k].m - exact_m(t));
    double a_error = fabs(history.states[k].a.i - exact_a(t));

    /*
    printf("%f,%f,%f,%f,", t, history.states[k].m, exact_m(t), m_error);
    printf("%f,%f,%f,", history.states[k].v.i, exact_v(t), v_error);
    printf("%f,%f,%f,", history.states[k].a.i, exact_a(t), a_error);
    printf("%f,%f,%f\n", altitude(history.states[k].x.i), exact_x(t), x_error);
    */
    
    char * err = "\n  (-) Error: thrust_test1()\n        (+) Integration error larger than expected\n";
    
    mu_assert(err, x_error < 0.001);
    mu_assert(err, v_error < 0.001);
    mu_assert(err, m_error < 0.00001);
    mu_assert(err, a_error < 0.001);
  }
  
  return 0;
}

/**
 * @test This sets up a 1D thrust and balistics model and compairs the result 
 * from astraea() to the true, analytical solution at all returned points.
 */
char * thrust_test2()
{ 
  init_rocket_setup();
  
  double intial_velocity = 0.0;
  double inital_height = 0.0;
  double Isp = 260;
  double m_0 = 25.4;
  double avg_thrust = 2500;
  double fuel = 4.5;
  
  // Exact
  double exact_x(double t)
  {
    double mdot = avg_thrust/(g_0*Isp);
    double mass = m_0 - mdot * t;
    
    double alt = g_0*Isp*(m_0/mdot)*(1-(mass/m_0)*(log(m_0/mass)+1)) - 0.5*g_0*t*t;
    
    //printf("m(x) = m - d*x;\n");
    //printf("f(x) = g*i*(m/d)*(1-(m(x)/m)*(log(m/m(x))+1)) - 0.5*g*x*x;\n");
    //printf("g = %f;i = %f;m = %f;d = %f;\n", g_0, Isp, m_0, mdot);
    
    return alt;
  }
  
  double exact_v(double t)
  {
    double mass = m_0 - avg_thrust/(g_0*Isp) * t;
    double velocity = g_0*Isp*log(m_0 / mass) - g_0*t + intial_velocity;
    return velocity;
  }
  
  double exact_m(double t)
  {
    double mass = m_0 - avg_thrust/(g_0*Isp) * t;
    return mass;
  }
  
  double exact_a(double t)
  {
    double mdot = avg_thrust/(g_0*Isp);
    double mass = m_0 - mdot * t;
    double acc = ((g_0*Isp)/mass) * mdot - g_0;
    return acc;
  }
  
  double exact_burn_time = fuel / (avg_thrust/(g_0*Isp));
  
  double exact_v_bo = exact_v(exact_burn_time);
  double exact_x_bo = exact_x(exact_burn_time);
  
  double exact_x_coast(double t)
  {
    t = t - exact_burn_time;
    double alt = -0.5*g_0*t*t + exact_v_bo*t + exact_x_bo;
    return alt;
  }
  
  //g(x) = -0.5*g*(x + t)*(x * t) + v*(x + t) + b;
  //printf("v = %f; b = %f;t = %f\n", exact_v_bo, exact_x_bo, exact_burn_time);
  
  double exact_v_coast(double t)
  {
    t = t - exact_burn_time;
    double vel = -g_0*t + exact_v_bo;
    return vel;
  }
  
  double sim_time;
  
  init.x.i += inital_height;
  init.v.i = intial_velocity;
  init.m = m_0;
  
  rocket_stage sustainer;
  sustainer.area = 0.1;
  sustainer.Cd = 0.8;
  
  thrust_curve th;
  build_thrust_curve(fuel, Isp, avg_thrust, &th);
  
  double burn_time = th.time[th.n];
 
  //printf("%f == %f\n\n",burn_time, exact_burn_time);
  char * err = "\n  (-) Error: thrust_test2()\n        (+) Incorrect Burn Time\n";
  double burn_time_error = fabs(burn_time - exact_burn_time);
  mu_assert(err, burn_time_error < 0.00001);
  
  sustainer.thrust = th;
  
  astraea(  init, 
            sustainer,
            0,
            burn_time,
            &history,
            &sim_time,
            physics_OneDOF_Fixed_G_Thrust_No_Atmosphere,
            boundary_condition_max_alt);
  
  
  astraea(  history.states[history.length - 1],
            sustainer,
            burn_time,
            2046,
            &history,
            &sim_time,
            physics_OneDOF_Fixed_G_Thrust_No_Atmosphere,
            boundary_condition_max_alt);
  
  
  err = "\n  (-) Error: thrust_test3()\n        (+) Integration error larger than expected\n";
  
  ///printf("%d\n",history.length);
  int k;
  for (k=0;k<history.length;k++)
  {
    double t =  history.times[k];
    
    if (t < exact_burn_time)
    {
      double x_error = fabs(altitude(history.states[k].x) - exact_x(t));
      double v_error = fabs(history.states[k].v.i - exact_v(t));
      double m_error = fabs(history.states[k].m - exact_m(t));
      double a_error = fabs(history.states[k].a.i - exact_a(t));
      
      /*
      printf("%f,%f,%f,%f", t, altitude(history.states[k].x), exact_x(t), x_error);
      printf("%f,%f,%f,", history.states[k].v.i, exact_v(t), v_error);
      printf("%f,%f,%f,", history.states[k].a.i, exact_a(t), a_error);
      
      printf("%f,%f,%f\n", history.states[k].m, exact_m(t), m_error);
      */
      
      mu_assert(err, x_error < 0.001);
      mu_assert(err, v_error < 0.001);
      mu_assert(err, m_error < 0.00001);
      mu_assert(err, a_error < 0.001);
    }
    
    if (t > exact_burn_time)
    {
      double x_error = fabs(altitude(history.states[k].x) - exact_x_coast(t));
      double v_error = fabs(history.states[k].v.i - exact_v_coast(t));
      double a_error = fabs(history.states[k].a.i - -g_0);
      
      /*
      printf("%f,%f,%f,%f,", t, altitude(history.states[k].x), exact_x_coast(t), x_error);
      printf("%f,%f,%f,", history.states[k].v.i, exact_v_coast(t), v_error);
      printf("%f,%f,%f,\n", history.states[k].a.i, -g_0, a_error);
      */
      
      mu_assert(err, x_error < 0.001);
      mu_assert(err, v_error < 0.001);
      mu_assert(err, a_error < 0.001);
    }
  }
  
  return 0;
}

void init_rocket_setup()
{
  /// Dummy values
  init.x.i = RADIUS_EARTH;
  init.x.j = 0.0;
  init.x.k = 0.0;
  init.v.i = 0.0;
  init.v.j = 0.0;
  init.v.k = 0.0;
  init.a.i = 0.0;
  init.a.j = 0.0;
  init.a.k = 0.0;
  init.m = 1.0;
  
  history.times = NULL;
  history.states = NULL;
  history.length = 0;
}
