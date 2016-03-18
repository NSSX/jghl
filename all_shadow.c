#include "my_h.h"

int shadow_sphere(t_sphere *sphere, t_light *light, t_ray *ray, t_struct *mystruct, double v)
{
  double a;
  double b;
  double c;
  double det;
  double myv1;
  double myv2;

  give_sphere(sphere, mystruct->lights, ray, v, mystruct);
  a = (mypow(mystruct->lights->dx, 2) + mypow(mystruct->lights->dy, 2) + mypow(mystruct->lights->dz, 2));
  b = 2 * (mystruct->lights->dx * mystruct->xp + mystruct->lights->dy * mystruct->yp
           + mystruct->lights->dz * mystruct->zp);
  c = (mypow(mystruct->xp, 2) + mypow(mystruct->yp, 2)
       + mypow(mystruct->zp, 2) - mypow(sphere->radius, 2));
  det = (b * b) - (4 * a * c);
  if(det < 0)
    return (1);
  myv1 = (-b - sqrt(det)) / (2 * a);
  myv2 = (-b + sqrt(det)) / (2 * a);
  if(myv1 > myv2 && myv2 > 0)
    {
      myv1 = myv2;
    }
  if(myv1 > 0 && myv1 < 1)
    return (2);
  return (1);
}

int shadow_cyl(t_sphere *cyl, t_light *light, double v, t_ray *ray)
{
  double a;
  double b;
  double c;
  double det;
  double myv1;
  double myv2;

  give_cylinder(cyl, light, ray, v);
  a = (mypow(light->dx, 2) + mypow(light->dy, 2));
  b = (2 * light->dx * cyl->xt + 2 * light->dy * cyl->yt);
  c = (mypow(cyl->xt, 2) + mypow(cyl->yt, 2) - cyl->radius * cyl->radius);
  det = (b * b) - (4 * a * c);
  if(det < 0)
    return (1);
  myv1 = (-b - sqrt(det)) / (2 * a);
  myv2 = (-b + sqrt(det)) / (2 * a);
  if(myv1 > myv2 && myv2 > 0)
    {
      myv1 = myv2;
    }
  if(myv1 > 0 && myv1 < 1)
    return (2);
  return (1);
}

int shadow_cone(t_sphere *cone, t_light *light, double v, t_ray *ray)
{
  double a;
  double b;
  double c;
  double det;
  double myv1;
  double myv2;

  give_cone(cone, light, ray, v);
  a = mypow(light->dx, 2) + mypow(light->dy, 2) - (mypow(light->dz, 2) * 0.00000001);
  b = (2 * light->dx * cone->xt) + (2 * light->dy * cone->yt) - ((2 * light->dz * cone->zt) * 0.00000001)\
    ;
  c = (mypow(cone->xt, 2) + mypow(cone->yt, 2) - (mypow(cone->zt, 2) * 0.00000001));
  det = (b * b) - (4 * a * c);
  if(det < 0)
    return (1);
  myv1 = (-b - sqrt(det)) / (2 * a);
  myv2 = (-b + sqrt(det)) / (2 * a);
  if(myv1 > myv2 && myv2 > 0)
    {
      myv1 = myv2;
    }
  if(myv1 > 0 && myv1 < 1)
    return (2);
  return (1);
}

