#include "my_h.h"

int             calc_shadow_sph(t_point *point, t_coord *sph, float k);
int             calc_shadow_cyl(t_point *point, t_coord *cyl, float k);
  int             calc_shadow_cone(t_point *point, t_coord *cone, float k);

int light_plan(t_point *point, t_struct *mystruct)
{
  double v1;
  double n1;
  double n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = point->coord->xeye + point->coord->vx * mystruct->plan->v;
  tempy = point->coord->yeye + point->coord->vy * mystruct->plan->v;
  tempz = point->coord->zeye + point->coord->vz * mystruct->plan->v;
  point->lum->vx = point->lum->x - tempx;
  point->lum->vy = point->lum->y - tempy;
  point->lum->vz = point->lum->z - tempz;
  v1 = (point->lum->nx * point->lum->vx + point->lum->ny * point->lum->vy + point->lum->nz * point->\
        lum->vz);
  n1 = (pow(point->lum->nx,2) + pow(point->lum->ny,2) + pow(point->lum->nz, 2));
  n2 = (pow(point->lum->vx,2) + pow(point->lum->vy,2) + pow(point->lum->vz, 2));
  mystruct->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  if (calc_shadow_sph(point, point->sph, mystruct->plan->v) == 2 ||
    calc_shadow_cyl(point, point->cyl, mystruct->plan->v) == 2 ||
     calc_shadow_cone(point, point->cone, mystruct->plan->v) == 2 )
    return (0x000000);
  return (0x00FF00);
}



void            my_shadow_sph(t_point *point, t_coord *sph, float k)
{
  point->pos->xp = point->coord->xeye + (sph->vx * k);
  point->pos->yp = point->coord->yeye + (sph->vy * k);
  point->pos->zp = point->coord->zeye + (sph->vz * k);
  point->lum->vx = point->lum->x - point->pos->xp;
  point->lum->vy = point->lum->y - point->pos->yp;
  point->lum->vz = point->lum->z - point->pos->zp;
}

void            my_shadow_cyl(t_point *point, t_coord *cyl, float k)
{
  point->pos->xp = point->coord->xeye + (cyl->vx * k);
  point->pos->yp = point->coord->yeye + (cyl->vy * k);
  point->lum->vx = point->lum->x - point->pos->xp;
  point->lum->vy = point->lum->y - point->pos->yp;
}

void            my_shadow_cone(t_point *point, t_coord *cone, float k)
{
  point->pos->xp = point->coord->xeye + (cone->vx * k);
  point->pos->yp = point->coord->yeye + (cone->vy * k);
  point->pos->zp = point->coord->zeye + (cone->vz * k);
  point->lum->vx = point->lum->x - point->pos->xp;
  point->lum->vy = point->lum->y - point->pos->yp;
  point->lum->vz = point->lum->z - point->pos->zp;
}


int             calc_shadow_cyl(t_point *point, t_coord *cyl, float k)
{
  double         a;
  double         b;
  double         c;
  double         det;
  double         myv1;
  double         myv2;

  my_shadow_cyl(point, cyl, k);
  a = (pow(point->lum->vx, 2) + pow(point->lum->vy, 2));
  b = (2 * point->lum->vx * point->pos->xp + 2
       * point->lum->vy * point->pos->yp);
  c = (pow(point->pos->xp, 2) + pow(point->pos->yp, 2)
       - cyl->rayon * cyl->rayon);
  det = (b * b) - (4 * a * c);
  myv1 = (-b - sqrt(det)) / (2 * a);
  myv2 = (-b + sqrt(det)) / (2 * a);
  if(myv1 > myv2 && myv2 > 0)
    myv1 = myv2;
  if(myv1 > 0 && myv1 < 1)
    return (2);
  return (1);
}

int             calc_shadow_cone(t_point *point, t_coord *cone, float k)
{
  double         a;
  double         b;
  double         c;
  double         det;
  double         myv1;
  double         myv2;

  my_shadow_cone(point, cone, k);
  a = pow(point->lum->vx, 2) + pow(point->lum->vy, 2)
    - (pow(point->lum->vz, 2) * 0.05);
  b = (2 * point->lum->vx * point->pos->xp) + (2 * point->lum->vy
					       * point->pos->yp)
    - ((2 * point->lum->vz * point->pos->zp) * 0.05);
  c = (pow(point->pos->xp, 2) + pow(point->pos->yp, 2)
       - (pow(point->pos->zp, 2) * 0.05));
  det = (b * b) - (4 * a * c);
  myv1 = (-b - sqrt(det)) / (2 * a);
  myv2 = (-b + sqrt(det)) / (2 * a);
  if(myv1 > myv2 && myv2 > 0)
    myv1 = myv2;
  if(myv1 > 0 && myv1 < 1)
    return (2);
  return (1);
}


int             my_light_cyl(t_point *point, t_struct *mystruct)
{
  float         v1;
  float         norme2;
  float         norme1;

  point->pos->xp = point->coord->xeye + (point->cyl->vx * mystruct->cyls->pos->v);
  point->pos->yp = point->coord->yeye + (point->cyl->vy * mystruct->cyls->pos->v);
  point->pos->zp = point->coord->zeye + (point->cyl->vz * mystruct->cyls->pos->v);
  point->lum->nx = point->pos->xp;
  point->lum->ny = point->pos->yp;
  point->lum->nz = 0;
  point->lum->vx = point->lum->x - point->pos->xp;
  point->lum->vy = point->lum->y - point->pos->yp;
  point->lum->vz = point->lum->z - point->pos->zp;
  v1 = (point->lum->nx * point->lum->vx + point->lum->ny
        * point->lum->vy + point->lum->nz * point->lum->vz);
  norme1 = (pow(point->lum->nx, 2) + pow(point->lum->ny, 2)
            + pow(point->lum->nz, 2));
  norme2 = (pow(point->lum->vx, 2) + pow(point->lum->vy, 2)
            + pow(point->lum->vz, 2));
  mystruct->img->cos = v1 / (sqrt(norme1) * sqrt(norme2));
  if (calc_shadow_sph(point, point->sph, mystruct->cyls->pos->v) == 2 ||
    calc_shadow_cone(point, point->cone, mystruct->cyls->pos->v) == 2)
    return (0x000000);
  return (0x0000FF);
}





int             my_light_cone(t_point *point, t_struct *mystruct)
{
  float         v1;
  float         norme2;
  float         norme1;

  point->pos->xp = point->coord->xeye + (point->cone->vx * mystruct->cones->pos->v);
  point->pos->yp = point->coord->yeye + (point->cone->vy * mystruct->cones->pos->v);
  point->pos->zp = point->coord->zeye + (point->cone->vz * mystruct->cones->pos->v);
  point->lum->nx = point->pos->xp;
  point->lum->ny = point->pos->yp;
  point->lum->nz = 0.05 * point->pos->zp;
  point->lum->vx = point->lum->x - point->pos->xp;
  point->lum->vy = point->lum->y - point->pos->yp;
  point->lum->vz = point->lum->z - point->pos->zp;
  v1 = (point->lum->nx * point->lum->vx + point->lum->ny
        * point->lum->vy + point->lum->nz * point->lum->vz);
  norme1 = (pow(point->lum->nx, 2) + pow(point->lum->ny, 2)
            + pow(point->lum->nz, 2));
  norme2 = (pow(point->lum->vx, 2) + pow(point->lum->vy, 2)
            + pow(point->lum->vz, 2));
  mystruct->img->cos = v1 / (sqrt(norme1) * sqrt(norme2));
  if (calc_shadow_sph(point, point->sph, mystruct->cones->pos->v) == 2 ||
      calc_shadow_cyl(point, point->cyl, mystruct->cones->pos->v) == 2)
    return (0x000000);
  return (0x0000FF);
}

int             my_light_sph(t_point *point, t_struct *mystruct)
{
  float         v1;
  float         norme2;
  float         norme1;
  double tempx;
  double tempy;
  double tempz;

  tempx = point->coord->xeye + (point->sph->vx * mystruct->spheres->pos->v);
  tempy = point->coord->yeye + (point->sph->vy * mystruct->spheres->pos->v);
  tempz = point->coord->zeye + (point->sph->vz * mystruct->spheres->pos->v);
  point->lum->nx = tempx;
  point->lum->ny = tempy;
  point->lum->nz = tempz;
  point->lum->vx = point->lum->x - tempx;
  point->lum->vy = point->lum->y - tempy;
  point->lum->vz = point->lum->z - tempz;
  v1 = (point->lum->nx * point->lum->vx + point->lum->ny
        * point->lum->vy + point->lum->nz * point->lum->vz);
  norme1 = (pow(point->lum->nx, 2) + pow(point->lum->ny, 2)
            + pow(point->lum->nz, 2));
  norme2 = (pow(point->lum->vx, 2) + pow(point->lum->vy, 2)
            + pow(point->lum->vz, 2));
  mystruct->img->cos = v1 / (sqrt(norme1) * sqrt(norme2));
  if (calc_shadow_cyl(point, point->cyl, mystruct->spheres->pos->v) == 2 ||
      calc_shadow_cone(point, point->cone, mystruct->spheres->pos->v) == 2)
    return (0x000000);
  return (0xFF0000);
}

