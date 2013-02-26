/**
 * norm
 */
double norm(vec v);

/**
 * unit_vec
 */
vec unit_vec(vec v);
vec vec_scale(vec v, double s);
double dot_prod(vec a, vec b);
vec matrix_mult(mat3 m, vec v);
mat3 axis_angle_to_rotation_matrix(vec axis_angle);
