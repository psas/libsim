/**
 * Main
 *
 * Author: Nathan Bergey
 * Version: Wednesday, April 13 2011
 */

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <time.h>
#include "utils/structs.h"
#include "libsim.h"
#include "physics/physics.h"
#include "physics/gravity.h"
#include "physics/models/earth.h"
#include "physics/models/thrustcurve.h"
#include "utils/coord.h"
#include "utils/boundary_conditions.h"


void switches(int argc, char **argv);
void print_history(void);
void print_line(FILE *out, double t, state s);
void build_thrust_curve(double fuel, double isp, double thrust, thrust_curve *t);
void OneDOF_Vertical_Sphereical_Earth_Gravity_Exp_Atmosphere_Height(double fuel, double isp, double avg_thrust, double begin_mass, double end_mass, double opt_resolution);
void OneDOF_Constant_Gravity_No_Atmosphere_Height(double fuel, double isp, double avg_thrust, double begin_mass, double end_mass, double opt_resolution);
void setup(void);
void OneDOF_vertical_single_run(void);
void OneDOF_vertical_thrust_single_run(double fuel);

float interp_dt = 0;
int do_print_history = 0;
char *outFileName = "out.dat";

/* Globals */
/// Initial State of the rocket.
state init;

/// State and time history of one sim run.
state_history history;

/// The Rocket
rocket the_rocket;

/// Simulation Run Time
double sim_time = 0;
double total_sim_time = 0;
double average_sim_time = 0;
double max_sim_time = 0;
double min_sim_time = 1000000;


/**
 * main
 */
int main(int argc, char **argv)
{
  /// Procces Switches
  switches(argc, argv);
  
  /// Initilize memory
  setup();
  
  OneDOF_vertical_single_run();
  //OneDOF_vertical_thrust_single_run(8.0);
  
  print_history();

  /// Exit
  return 0;
}

void switches(int argc, char **argv)
{
  int i;
  /* Start at i = 1 to skip the command name. */
  for (i = 1; i < argc; i++) 
  {
    /* Check for a switch (leading "-"). */
    if (argv[i][0] == '-') 
    {
        /* Use the next character to decide what to do. */
        switch (argv[i][1]) 
        {
	        case 'i':   // interpolate results
	          interp_dt = atof(argv[i+1]);
			      break;
			    case 'h':   // print history
			      do_print_history = 1;
			      break;
			    case 'f':   // output file name
			      outFileName = argv[i+1];
			      break;
	        default:	
	          printf("Unknown switch %s\n", argv[i]);
	          exit(1);
        }
    }
  }
  return;
}

void print_history(void)
{
  int i;
  
  FILE *out;
  out = fopen(outFileName, "w");
  if ( out == NULL )
  {
      printf("File Handle error.\n");
      exit(1);
  }
  
  for (i=0;i<history.length;i++)
  {
    print_line(out, history.times[i], history.states[i]);
  }

  /*
    int k;
    FILE *out;
    out = fopen(outFileName, "w");
    if ( out == NULL )
    {
        printf("File Handle error.\n");
        exit(1);
    }

    if (interp_dt > 0)
    {
      state history_i[MAXSTEPS];
      interpolate_history(time, history, steps, interp_dt, history_i);
      double x = time[0];
      double xf = time[steps-2];
      int steps_to_take = ceil((xf-x)/interp_dt);
      for (k=0;k<=steps_to_take;k++)
      {
        print_line(out, x, history_i[k]);
        if ( (x + interp_dt) > xf)
          x = xf;
        else
          x += interp_dt;
      }
    }
    else
    {
      for (k=0;k<steps;k++)
      {
        print_line(out, time[k], history[k]);
      }
    }
    */
}

void print_line(FILE *out, double t, state s)
{
  char format[512] = "";
  char exp[8] = "%0.10e,";
  
  strcat(format, "%f,");
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, exp);
  strcat(format, "\n");
  
  fprintf(out, format
    , t                                     //  1
    , s.x.v.i                                 //  2
    , s.x.v.j                                 //  3
    , s.x.v.k                                 //  4
    , altitude(s.x)                         //  5
    , s.v.v.i                                 //  6
    , s.v.v.j                                 //  7
    , s.v.v.k                                 //  8
    , vertical_velocity(s)                  //  9
    , s.a.v.i                                 // 10
    , s.a.v.j                                 // 11
    , s.a.v.k                                 // 12
    , vertical_acceleration_gee(s.a)        // 13
    , s.m);                                 // 14
    
}

void setup(void)
{
  /// Dummy values
  init.x.v.i = RADIUS_EARTH;
  init.x.v.j = 0.0;
  init.x.v.k = 0.0;
  init.v.v.i = 0.0;
  init.v.v.j = 0.0;
  init.v.v.k = 0.0;
  init.a.v.i = 0.0;
  init.a.v.j = 0.0;
  init.a.v.k = 0.0;
  init.m = 1.0;
  
  history.length = 0;
  history.times = NULL;
  history.states = NULL;
}

void OneDOF_vertical_single_run(void)
{
  /// Setup
  init.v.v.i = 2390.0;
  init.m = 25.5;
  
  rocket_stage sustainer;
  
  integration_strategy strategy;
  
  strategy.DOF            = 3;
  strategy.boundary       = boundary_condition_max_alt;
  strategy.gravity_model  = gravity_sphere;
  
  printf("starting\n");
  
  /// Run
  history = astraea(init, 
                    sustainer,
                    0, 
                    1000,
                    &sim_time,
                    strategy);
  
  int i;
  for (i=0;i<history.length;i++) 
  {
    printf("%f,", history.times[i]);
    printf("%f\n", history.states[i].x.v.i);
  }
  
  /// Print status
  printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("Number of Steps:  %d\n", history.length);
  printf("Sim Time:         %0.3f ms\n", sim_time/1000.0);
  return;
}


void OneDOF_vertical_thrust_single_run(double fuel)
{
 /*
  /// Setup
  init.v.i = 0.0;
  init.m = 25.5 + fuel;
  
  rocket_stage sustainer;
  sustainer.area = 0.1;
  sustainer.Cd = 0.8;
  
  thrust_curve th;
  build_thrust_curve(fuel, 245, 3800, &th);
  
  double burn_time = th.time[th.n];
  
  sustainer.thrust = th;
  
  /// Run
  astraea(  init,
            sustainer,
            0, 
            burn_time,
            &history,
            &sim_time,
            physics_OneDOF_Fixed_G_Thrust_No_Atmosphere,
            boundary_condition_max_alt);
  
  
  //free(sustainer.thrust.time);
  //free(sustainer.thrust.m_dot);
  //sustainer.thrust.time = NULL;
  //sustainer.thrust.m_dot = NULL;
  //sustainer.thrust.n = 0;
  //sustainer = NULL;
  
  //printf("%f,%f\n", history.times[history.length], history.times[history.length - 1]);
  
  
  astraea(  history.states[history.length - 1],
            sustainer,
            burn_time,
            2046,
            &history,
            &sim_time,
            physics_OneDOF_Fixed_G_No_Atmosphere,
            boundary_condition_max_alt);
  
  */
  
  /// Print status
  printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("Number of Steps:  %d\n", history.length);
  printf("Sim Time:         %0.0f μs\n", sim_time);
  return;
}
