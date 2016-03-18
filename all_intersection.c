#include "my_h.h"

//int plan_int(t_ray *ray, t_vec3d *plan, t_light *light)
int plan_int(t_struct *mystruct)
{
 mystruct->plan->v = (-mystruct->rays->o->z) / mystruct->rays->d->z;
  if(mystruct->rays->d->z == 0)
    return (0);
  mystruct->lights->ox = 0;
  mystruct->lights->oy = 0;
  mystruct->lights->oz = 100;
  if(mystruct->plan->v > 0)
    return (2);
  else
    return (1);
}

void sphere_int(t_sphere *sphere, t_ray *ray)
{
  double a;
  double b;
  double c;
  double det;
  double myv;

  a = (sphere->pos->x * sphere->pos->x + sphere->pos->y * sphere->pos->y + sphere->pos->z *sphere->pos->z);
  b = (2 * sphere->pos->x * ray->o->x + 2 * sphere->pos->y * ray->o->y + 2 * sphere->pos->z * ray->o->z);
  c = (ray->o->x * ray->o->x + ray->o->y * ray->o->y + ray->o->z * ray->o->z - sphere->radius * sphere->radius);
  det = (b * b) - (4 * a * c);
  sphere->pos->v = (-b - sqrt(det)) / (2 * a);
  myv = (-b + sqrt(det)) / (2 * a);
  if(det < 0)
    sphere->pos->v = -1;
  else if (sphere->pos->v > myv && myv > 0)
    sphere->pos->v = myv;
}

void cone_int(t_sphere *cone, t_ray *ray)
{
  double a;
  double b;
  double c;
  double det;
  double myv;

  a = mypow(cone->pos->x, 2) + mypow(cone->pos->y, 2) - (mypow(cone->pos->z, 2) * 0.1);
  b = (2 * cone->pos->x * ray->o->x) + (2 * cone->pos->y * ray->o->y) - ((2 * cone->pos->z * ray->o->z) * 0.1);
  c = (mypow(ray->o->x, 2) + mypow(ray->o->y, 2) - (mypow(ray->o->z, 2) * 0.1));
  det = (b * b) - (4 * a * c);
  cone->pos->v = (-b - sqrt(det)) / (2 * a);
  myv = (-b + sqrt(det)) / (2 * a);
  if(det < 0)
    cone->pos->v = -1;
  else if (cone->pos->v > myv && myv > 0)
    cone->pos->v = myv;
}

void cyl_int(t_sphere *cyl, t_ray *ray)
{
  double a;
  double b;
  double c;
  double det;
  double myv;

  a = mypow(cyl->pos->x, 2) + mypow(cyl->pos->y, 2);
  b = (2 * cyl->pos->x * ray->o->x + 2 * cyl->pos->y * ray->o->y);
  c = (mypow(ray->o->x, 2) + mypow(ray->o->y, 2) - mypow(cyl->radius, 2));
  det = (b * b) - (4 * a * c);
  cyl->pos->v = (-b - sqrt(det)) / (2 * a);
  myv = (-b + sqrt(det)) / (2 * a);
  if(det < 0)
    cyl->pos->v = -1;
  else if (cyl->pos->v > myv && myv > 0)
    cyl->pos->v = myv;
}

