vec ECEF2GEO(vec v);
vec GEO2ECEF(vec v);
vec ECEF2ENU(vec v, double lon, double lat);
vec ENU2ECEF(vec enu, double lon, double lat);
double altitude(vec ecef);
double vertical_velocity(state r);
double vertical_acceleration_gee(vec a);
