#include "my_h.h"

void vector_normalize(t_vec3d *v)
{
  double id;
  id =  1/sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
  v->x = v->x * id;
  v->y = v->y *id;
  v->z = v->z *id;
}

double vector_dot(t_vec3d *a, t_vec3d *b)
{
  double res;

  res = (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
  return (res);
}

t_vec3d *vector_copy(t_vec3d *a)//bad type?                                                               
{
  t_vec3d *v;
  v = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  v->x = a->x;
  v->y = a->y;
  v->z = a->z;
  return (v);
}

t_vec3d *vector_sub(t_vec3d *a, t_vec3d *b)
{
  t_vec3d *v;
  v = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  v->x = a->x - b->x;
  v->y = a->y - b->y;
  v->z = a->z - b->z;
  return (v);
}

t_vec3d *vector_mul(t_vec3d *a, t_vec3d *b)
{
  t_vec3d *v;
  v = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  v->x = a->x * b->x;
  v->y = a->y * b->y;
  v->z = a->z * b->z;
  return (v);
}
