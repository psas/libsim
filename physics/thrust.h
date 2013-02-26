vec thrust(state s, double t, double *mdot);
void set_thrust_curve(thrust_curve thrust);
void build_thrust_curve(double fuel, double isp, double avg_thrust, thrust_curve *t);
