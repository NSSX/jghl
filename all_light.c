#include "my_h.h"

int light_cyl(t_sphere *cyl, t_light *light, t_ray *ray, t_struct *mystruct)
{
  double v1;
  double n1;
  double n2;
  double tempx;
  double tempy;
  double tempz;

  cyl->xt = ray->o->x + (cyl->pos->x * cyl->pos->v);
  cyl->yt = ray->o->y + (cyl->pos->y * cyl->pos->v);
  cyl->zt = ray->o->z + (cyl->pos->z * cyl->pos->v);
  light->ox = cyl->xt;
  light->oy = cyl->yt;
  light->oz = 0;
  light->dx = light->x - cyl->xt;
  light->dy = light->y - cyl->yt;
  light->dz = light->z - cyl->zt;
  v1 = (light->ox * light->dx + light->oy * light->dy + light->oz * light->dz);
  n1 = (mypow(light->ox,2) + mypow(light->oy,2) + mypow(light->oz, 2));
  n2 = (mypow(light->dx,2) + mypow(light->dy,2) + mypow(light->dz, 2));
  mystruct->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  // if(shadow_sphere(mystruct->spheres, light, cyl->pos->v, ray) == 2 ||                                   
  //shadow_cone(mystruct->cones, light, cyl->pos->v, ray) == 2)                                           
  // return (0x000000);                                                                                    
  return (0x0000FF);
}

int light_cone(t_sphere *cone, t_light *light, t_ray *ray, t_struct *mystruct)
{
  double v1;
  double n1;
  double n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = ray->o->x + (cone->pos->x * cone->pos->v);
  tempy = ray->o->y + (cone->pos->y * cone->pos->v);
  tempz = ray->o->z + (cone->pos->z * cone->pos->v);
  light->ox = tempx;
  light->oy = tempy;
  light->oz = 0.00000001 * tempz;
  light->dx = light->x - tempx;
  light->dy = light->y - tempy;
  light->dz = light->z - tempz;
  v1 = (light->ox * light->dx + light->oy * light->dy + light->oz * light->dz);
  n1 = (mypow(light->ox,2) + mypow(light->oy,2) + mypow(light->oz, 2));
  n2 = (mypow(light->dx,2) + mypow(light->dy,2) + mypow(light->dz, 2));
  mystruct->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  //if(shadow_sphere(mystruct->spheres, light, cone->pos->v, ray) == 2 ||                                 
  // shadow_cyl(mystruct->cyls, light, cone->pos->v, ray) == 2)                                           
  // return (0x000000);                                                                                   
  return (0x0000FF);
}

int light_sphere(t_sphere *sphere, t_light *light, t_ray *ray, t_struct *mystruct)
{
  double v1;
  double n1;
  double n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = ray->o->x + (sphere->pos->x * sphere->pos->v);
  tempy = ray->o->y + (sphere->pos->y * sphere->pos->v);
  tempz = ray->o->z + (sphere->pos->z * sphere->pos->v);
  light->ox = tempx;
  light->oy = tempy;
  light->oz = tempz;
  light->dx = light->x - tempx;
  light->dy = light->y - tempy;
  light->dz = light->z - tempz;
  v1 = (light->ox * light->dx + light->oy * light->dy + light->oz * light->dz);
  n1 = (mypow(light->ox,2) + mypow(light->oy,2) + mypow(light->oz, 2));
  n2 = (mypow(light->dx,2) + mypow(light->dy,2) + mypow(light->dz, 2));
  mystruct->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  // if(shadow_cyl(mystruct->cyls, light, sphere->pos->v, ray) == 2 ||                                    
  //  shadow_cone(mystruct->cones, light, sphere->pos->v, ray) == 2)                                      
  // return (0x000000);                                                                                   
  return (0xFF0000);
}
