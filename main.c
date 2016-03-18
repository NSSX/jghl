#include "my_h.h"

void            my_shadow_sph(t_point *point, t_coord *sph, float k);
int             my_light_sph(t_point *point, t_struct *mystruct);
int             my_light_cone(t_point *point, t_struct *mystruct);
int             my_light_cyl(t_point *point, t_struct *mystruct);
void            my_value(t_point *point, int x, int y);

int             calc_shadow_sph(t_point *point, t_coord *sph, float k)
{
  double a;
  double b;
  double c;
  double det;
  double myv1;
  double myv2;

  my_shadow_sph(point, sph, k);
  a = (pow(point->lum->vx, 2) + pow(point->lum->vy, 2) + pow(point->lum->vz, 2));
  b = 2 * (point->lum->vx * point->pos->xp + point->lum->vy * point->pos->yp
           + point->lum->vz * point->pos->zp);
  c = (pow(point->pos->xp, 2) + pow(point->pos->yp, 2)
       + pow(point->pos->zp, 2) - pow(sph->rayon, 2));
  det = (b * b) - (4 * a * c);
  if(det < 0)
    return (1);
  myv1 = (-b - sqrt(det)) / (2 * a);
  myv2 = (-b + sqrt(det)) / (2 * a);
  if(myv1 > myv2 && myv2 > 0)
      myv1 = myv2;
  if(myv1 > 0 && myv1 < 1)
    return (2);
  return (1);
}


int plan_inter(t_main *main)
{
  main->plan->v = (-main->eye->zeye) / main->eye->vz;
  if(main->eye->vz == 0)
    return (0);
  main->lum->nx = 0;
  main->lum->ny = 0;
  main->lum->nz = 100;
  if(main->plan->v > 0)
    return (2);
  else
    return (1);
}

void sphere_inter(t_main *main)
{
  double a;
  double b;
  double c;
  double det;
  double myv;

  a = (main->sphere->x * main->sphere->x + main->sphere->y * main->sphere->y + main->sphere->z * main->sphere->z);
  b = (2 * main->sphere->x * main->eye->xeye + 2 * main->sphere->y * main->eye->yeye + 2 * main->sphere->z * main->eye->zeye);
  c = (main->eye->xeye * main->eye->xeye + main->eye->yeye * main->eye->yeye + main->eye->zeye * main->eye->zeye - main->sphere->rayon * main->sphere->rayon);
  det = (b * b) - (4 * a * c);
  main->sphere->v = (-b - sqrt(det)) / (2 * a);
  myv = (-b + sqrt(det)) / (2 * a);
  if(det < 0)
    main->sphere->v = -1;
  else if (main->sphere->v > myv && myv > 0)
    main->sphere->v = myv;
}

void cone_inter(t_main *main)
{
  double a;
  double b;
  double c;
  double det;
  double myv;

  a = mypow(main->cone->x, 2) + mypow(main->cone->y, 2) - (mypow(main->cone->z, 2) * 0.1);
  b = (2 * main->cone->x * main->eye->xeye) + (2 * main->cone->y * main->eye->yeye) - ((2 * main->cone->z * main->eye->zeye) * 0.1);
  c = (mypow(main->eye->xeye, 2) + mypow(main->eye->yeye, 2) - (mypow(main->eye->zeye, 2) * 0.1));
  det = (b * b) - (4 * a * c);
  main->cone->v = (-b - sqrt(det)) / (2 * a);
  myv = (-b + sqrt(det)) / (2 * a);
  if(det < 0)
    main->cone->v = -1;
  else if (main->cone->v > myv && myv > 0)
    main->cone->v = myv;
}

void cyl_inter(t_main *main)
{
  double a;
  double b;
  double c;
  double det;
  double myv;

  a = mypow(main->cyl->x, 2) + mypow(main->cyl->y, 2);
  b = (2 * main->cyl->x * main->eye->xeye + 2 * main->cyl->y * main->eye->yeye);
  c = (mypow(main->eye->xeye, 2) + mypow(main->eye->yeye, 2) - mypow(main->cyl->rayon, 2));
  det = (b * b) - (4 * a * c);
  main->cyl->v = (-b - sqrt(det)) / (2 * a);
  myv = (-b + sqrt(det)) / (2 * a);
  if(det < 0)
    main->cyl->v = -1;
  else if (main->cyl->v > myv && myv > 0)
    main->cyl->v = myv;
}

void            give_shadow_sphere(t_main *main, double v)
{
  main->sphere->xp = main->eye->xeye + (main->sphere->x * v);
  main->sphere->yp = main->eye->yeye + (main->sphere->y * v);
  main->sphere->zp = main->eye->zeye + (main->sphere->z * v);
  main->lum->vx = main->lum->x - main->sphere->xp;
  main->lum->vy = main->lum->y - main->sphere->yp;
  main->lum->vz = main->lum->z - main->sphere->zp;
}

void            give_shadow_cyl(t_main *main, double v)
{
  main->cyl->xp = main->eye->xeye + (main->cyl->x * v);
  main->cyl->yp = main->eye->yeye + (main->cyl->y * v);
  main->lum->vx = main->lum->x - main->cyl->xp;
  main->lum->vy = main->lum->y - main->cyl->yp;
}

void            give_shadow_cone(t_main *main, double v)
{
  main->cone->xp = main->eye->xeye + (main->cone->x * v);
  main->cone->yp = main->eye->yeye + (main->cone->y * v);
  main->cone->zp = main->eye->zeye + (main->cone->z * v);
  main->lum->vx = main->lum->x - main->cone->xp;
  main->lum->vy = main->lum->y - main->cone->yp;
  main->lum->vz = main->lum->z - main->cone->zp;
}

int             my_shadow_cylinder(t_main *main, double v)
{
  double         a;
  double         b;
  double         c;
  double         det;
  double         myv1;
  double         myv2;

  give_shadow_cyl(main, v);
  a = (pow(main->lum->vx, 2) + pow(main->lum->vy, 2));
  b = (2 * main->lum->vx * main->cyl->xp + 2
       * main->lum->vy * main->cyl->yp);
  c = (pow(main->cyl->xp, 2) + pow(main->cyl->yp, 2)
       - main->cyl->rayon * main->cyl->rayon);
  det = (b * b) - (4 * a * c);
  myv1 = (-b - sqrt(det)) / (2 * a);
  myv2 = (-b + sqrt(det)) / (2 * a);
  if(myv1 > myv2 && myv2 > 0)
    myv1 = myv2;
  if(myv1 > 0 && myv1 < 1)
    return (2);
  return (1);
}

int             my_shadow_conee(t_main *main, double v)
{
  double         a;
  double         b;
  double         c;
  double         det;
  double         myv1;
  double         myv2;

  give_shadow_cone(main, v);
  a = pow(main->lum->vx, 2) + pow(main->lum->vy, 2)
    - (pow(main->lum->vz, 2) * 0.05);
  b = (2 * main->lum->vx * main->cone->xp) + (2 * main->lum->vy
                                               * main->cone->yp)
    - ((2 * main->lum->vz * main->cone->zp) * 0.05);
  c = (pow(main->cone->xp, 2) + pow(main->cone->yp, 2)
       - (pow(main->cone->zp, 2) * 0.05));
  det = (b * b) - (4 * a * c);
  myv1 = (-b - sqrt(det)) / (2 * a);
  myv2 = (-b + sqrt(det)) / (2 * a);
  if(myv1 > myv2 && myv2 > 0)
    myv1 = myv2;
  if(myv1 > 0 && myv1 < 1)
    return (2);
  return (1);
}

int             my_shadow_sphere(t_main *main, double v)
{
  double a;
  double b;
  double c;
  double det;
  double myv1;
  double myv2;

  give_shadow_sphere(main, v);
  a = (pow(main->lum->vx, 2) + pow(main->lum->vy, 2) + pow(main->lum->vz, 2));
  b = 2 * (main->lum->vx * main->sphere->xp + main->lum->vy * main->sphere->yp
           + main->lum->vz * main->sphere->zp);
  c = (pow(main->sphere->xp, 2) + pow(main->sphere->yp, 2)
       + pow(main->sphere->zp, 2) - pow(main->sphere->rayon, 2));
  det = (b * b) - (4 * a * c);
  if(det < 0)
    return (1);
  myv1 = (-b - sqrt(det)) / (2 * a);
  myv2 = (-b + sqrt(det)) / (2 * a);
  if(myv1 > myv2 && myv2 > 0)
    myv1 = myv2;
  if(myv1 > 0 && myv1 < 1)
    return (2);
  return (1);
}

int lighte_plan(t_point *point, t_main *main)
{
  double v1;
  double n1;
  double n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = point->coord->xeye + point->coord->vx * main->plan->v;
  tempy = point->coord->yeye + point->coord->vy * main->plan->v;
  tempz = point->coord->zeye + point->coord->vz * main->plan->v;
  point->lum->vx = point->lum->x - tempx;
  point->lum->vy = point->lum->y - tempy;
  point->lum->vz = point->lum->z - tempz;
  v1 = (point->lum->nx * point->lum->vx + point->lum->ny * point->lum->vy + point->lum->nz * point->\
        lum->vz);
  n1 = (pow(point->lum->nx,2) + pow(point->lum->ny,2) + pow(point->lum->nz, 2));
  n2 = (pow(point->lum->vx,2) + pow(point->lum->vy,2) + pow(point->lum->vz, 2));
  main->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  if (calc_shadow_sph(point, point->sph, main->plan->v) == 2)
    return (0x000000);
  return (0x00FF00);
}

int my_light_plan(t_main *main)
{
  double v1;
  double n1;
  double n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = main->eye->xeye + main->eye->vx * main->plan->v;
  tempy = main->eye->yeye + main->eye->vy * main->plan->v;
  tempz = main->eye->zeye + main->eye->vz * main->plan->v;
  main->lum->vx = main->lum->x - tempx;
  main->lum->vy = main->lum->y - tempy;
  main->lum->vz = main->lum->z - tempz;
  v1 = (main->lum->nx * main->lum->vx + main->lum->ny * main->lum->vy + main->lum->nz * main->\
        lum->vz);
  n1 = (pow(main->lum->nx,2) + pow(main->lum->ny,2) + pow(main->lum->nz, 2));
  n2 = (pow(main->lum->vx,2) + pow(main->lum->vy,2) + pow(main->lum->vz, 2));
  main->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  if (my_shadow_sphere(main, main->plan->v) == 2 ||
    my_shadow_cylinder(main, main->plan->v) == 2 ||
      my_shadow_conee(main, main->plan->v) == 2 )
    return (0x000000);
  return (0x00FFFF);
}



int             my_light_cylinder(t_main *main)
{
  double v1;
  double n1;
  double n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = main->eye->xeye + (main->cyl->x * main->cyl->v);
  tempy = main->eye->yeye + (main->cyl->y * main->cyl->v);
  tempz = main->eye->zeye + (main->cyl->z * main->cyl->v);
  main->lum->nx = tempx;
  main->lum->ny = tempy;
  main->lum->nz = 0;
  main->lum->vx = main->lum->x - tempx;
  main->lum->vy = main->lum->y - tempy;
  main->lum->vz = main->lum->z - tempz;
  v1 = (main->lum->nx * main->lum->vx + main->lum->ny
        * main->lum->vy + main->lum->nz * main->lum->vz);
  n1 = (pow(main->lum->nx, 2) + pow(main->lum->ny, 2)
            + pow(main->lum->nz, 2));
  n2 = (pow(main->lum->vx, 2) + pow(main->lum->vy, 2)
            + pow(main->lum->vz, 2));
  main->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  if (my_shadow_sphere(main, main->cyl->v) == 2 ||
      my_shadow_conee(main, main->cyl->v) == 2)
	return (0x000000);
    return (0x0000FF);
}

int             my_light_conee(t_main *main)
{
  float         v1;
  float         n1;
  float         n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = main->eye->xeye + (main->cone->x * main->cone->v);
  tempy = main->eye->yeye + (main->cone->y * main->cone->v);
  tempz = main->eye->zeye + (main->cone->z * main->cone->v);
  main->lum->nx = tempx;
  main->lum->ny = tempy;
  main->lum->nz = 0.05 * tempz;
  main->lum->vx = main->lum->x - tempx;
  main->lum->vy = main->lum->y - tempy;
  main->lum->vz = main->lum->z - tempz;
  v1 = (main->lum->nx * main->lum->vx + main->lum->ny
        * main->lum->vy + main->lum->nz * main->lum->vz);
  n1 = (pow(main->lum->nx, 2) + pow(main->lum->ny, 2)
            + pow(main->lum->nz, 2));
  n2 = (pow(main->lum->vx, 2) + pow(main->lum->vy, 2)
            + pow(main->lum->vz, 2));
  main->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  if (my_shadow_sphere(main, main->cone->v) == 2 ||
      my_shadow_cylinder(main, main->cone->v) == 2)
      return (0x000000);
  return (0x0000FF);
}

int             my_light_sphere(t_main *main)
{
  float         v1;
  float         n1;
  float         n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = main->eye->xeye + (main->sphere->x * main->sphere->v);
  tempy = main->eye->yeye + (main->sphere->y * main->sphere->v);
  tempz = main->eye->zeye + (main->sphere->z * main->sphere->v);
  main->lum->nx = tempx;
  main->lum->ny = tempy;
  main->lum->nz = tempz;
  main->lum->vx = main->lum->x - tempx;
  main->lum->vy = main->lum->y - tempy;
  main->lum->vz = main->lum->z - tempz;
  v1 = (main->lum->nx * main->lum->vx + main->lum->ny
        * main->lum->vy + main->lum->nz * main->lum->vz);
  n1 = (pow(main->lum->nx, 2) + pow(main->lum->ny, 2)
            + pow(main->lum->nz, 2));
  n2 = (pow(main->lum->vx, 2) + pow(main->lum->vy, 2)
            + pow(main->lum->vz, 2));
  main->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  if (my_shadow_cylinder(main, main->sphere->v) == 2 ||
      my_shadow_conee(main, main->sphere->v) == 2)
     return (0x000000);
  return (0xFF0000);
}

void define_main(int x, int y, t_main *main)
{
  main->lum->nx = 0;
  main->lum->ny = 0;
  main->lum->nz = 100;
  main->lum->x = -3000;
  main->lum->y = 0;
  main->lum->z = 2000;
  main->eye->vx = 0;
  main->eye->vy = (WIDTH / 2) - x;
  main->eye->vz = (HEIGHT / 2) - y;
  main->eye->xeye = -3000;
  main->eye->yeye = 0;
  main->eye->zeye = 500;
  main->sphere->x = main->eye->vx + 1000- 400;
  main->sphere->y = main->eye->vy + 0;
  main->sphere->z = main->eye->vz + -300;
  main->sphere->rayon = 800;
  main->cyl->x = main->eye->vx + 1000;
  main->cyl->y = main->eye->vy + 0 - 100;
  main->cyl->z = main->eye->vz + -300;
  main->cyl->rayon = 150;
  main->cone->x = main->eye->vx + 1000;
  main->cone->y = main->eye->vy + 100;
  main->cone->z = main->eye->vz + -300;
}

void    my_give_v(double *v_value, t_main *main)
{
  v_value[0] = main->plan->v;
  v_value[1] = main->cone->v;
  v_value[2] = main->cyl->v;
  v_value[3] = main->sphere->v;
}

int color_choice(t_main *main)
{
  double v;
  double v_value[4];
  int i;
  double temp;


  v = 0;
  i = -1;
  my_give_v(v_value, main);
  while(i < 2)
    {
      i++;
      if(v_value[i] > v_value[i + 1])
        {
          temp = v_value[i];
          v_value[i] = v_value[i + 1];
          v_value[i + 1] = temp;
          i = -1;
        }
    }
  i = 0;
  while(i < 4 && v_value[i] <= 0)
    {
      i++;
      if(i == 4)
        break;
    }
  if(v > 3)
    v = v_value[3];
  else
    v = v_value[i];
  my_give_v(v_value, main);
  if(v == v_value[0])
    return (my_light_plan(main));
   else if(v == v_value[1])
    return(my_light_conee(main));
  else if(v == v_value[2])
    return(my_light_cylinder(main));
  else if(v == v_value[3])
    return(my_light_sphere(main));
  return (0x000000);
}



int definee_color(int x, int y, t_main *main)
{
   int color;
   
   main->sphere = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  main->cyl = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  main->cone = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  main->lum = (t_lum *)malloc(sizeof(t_lum) * 1);
  main->eye = (t_coord *)malloc(sizeof(t_coord) * 1);
  main->plan = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  define_main(x, y, main);
  plan_inter(main);
  sphere_inter(main);
  cone_inter(main);
  cyl_inter(main);
  color = color_choice(main);//start the lets choice func;
  return (color);
}



void            my_value(t_point *point, int x, int y)
{
  point->lum->nx = 0;
  point->lum->ny = 0;
  point->lum->nz = 100;
  point->lum->x = -3000;
  point->lum->y = 0;
  point->lum->z = 2000;
  point->coord->vx = 0;
  point->coord->vy = (WIDTH / 2) - x;
  point->coord->vz = (HEIGHT / 2) - y;
  point->coord->xeye = -3000;
  point->coord->yeye = 0;
  point->coord->zeye = 500;
  point->sph->vx = point->coord->vx + 1000- 400;
  point->sph->vy = point->coord->vy + 0;
  point->sph->vz = point->coord->vz + -300;
  point->sph->rayon = 800;
  point->cyl->vx = point->coord->vx + 1000;
  point->cyl->vy = point->coord->vy + 0 - 100;
  point->cyl->vz = point->coord->vz + -300;
  point->cyl->rayon = 150;
  point->cone->vx = point->coord->vx + 1000;
  point->cone->vy = point->coord->vy + 100;
  point->cone->vz = point->coord->vz + -300;
}

int lets_choice(t_sphere *sphere, t_vec3d *plan, t_light *light, t_ray *ray, t_struct *mystruct, t_sphere *cone, t_sphere *cyl, t_point *point)
{
   double v;
  double v_value[4];
  int i;
  double temp;


  v = 0;
  i = -1;
  give_v(plan, sphere, cone, cyl, v_value, point);
  while(i < 2)
    {
      i++;
      if(v_value[i] > v_value[i + 1])
        {
          temp = v_value[i];
          v_value[i] = v_value[i + 1];
          v_value[i + 1] = temp;
          i = -1;
        }
    }
  i = 0;
  while(i < 4 && v_value[i] <= 0)
    {
      i++;
      if(i == 4)
	break;
    }
  if(v > 3)
    v = v_value[3];
  else
    v = v_value[i];
  give_v(mystruct->plan, sphere, cone, cyl, v_value, point);
  if(v == v_value[0])
    return (light_plan(point, mystruct));
  else if(v == v_value[1])
       return(my_light_cone(point, mystruct));
  else if(v == v_value[2])
    return(my_light_cyl(point, mystruct));
  else if(v == v_value[3])
    return(my_light_sph(point, mystruct));
  return (0x000000);
}



int define_color(int x, int y, t_struct *mystruct)
{
  int color;
  t_sphere *sphere;
  t_vec3d *plan;
  t_ray *ray;
  t_light *light;
  t_sphere *cone;
  t_sphere *cyl;
  t_point *point;

  point = (t_point *)malloc(sizeof(t_point) * 1);
  light = (t_light *)malloc(sizeof(t_light) * 1);
  ray = (t_ray *)malloc(sizeof(t_ray) * 1);
  ray->o = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  ray->d = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  plan = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  cone = (t_sphere *)malloc(sizeof(t_sphere) * 1);
  cone->pos = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  sphere = (t_sphere *)malloc(sizeof(t_sphere) * 1);
  sphere->pos = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  cyl = (t_sphere *)malloc(sizeof(t_sphere) * 1);
  cyl->pos = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  point->img = (t_img *)malloc(sizeof(t_img) * 1);
  point->coord = (t_coord *)malloc(sizeof(t_coord) * 1);
  point->sph = (t_coord *)malloc(sizeof(t_coord) * 1);
  point->cyl = (t_coord *)malloc(sizeof(t_coord) * 1);
  point->cone = (t_coord *)malloc(sizeof(t_coord) * 1);
  point->pos = (t_pos *)malloc(sizeof(t_pos) * 1);
  point->lum = (t_lum *)malloc(sizeof(t_lum) * 1);
  mystruct->plan = plan;
  mystruct->rays = ray;
  mystruct->spheres = sphere;
  mystruct->cones = cone;
  mystruct->cyls = cyl;
  mystruct->lights = light;
  mystruct->point = point;
  point->count = mystruct->count;
  color = 0;
  give_value(x, y,sphere, ray, mystruct->lights, mystruct, cone, cyl);
  my_value(point,x,y);
  mystruct->rays->o->y = point->coord->yeye;
  plan_int(mystruct);  
  sphere_int(mystruct->spheres, ray);
  cone_int(mystruct->cones, ray);
  cyl_int(mystruct->cyls,ray);
  color = lets_choice(sphere, plan,mystruct->lights, ray, mystruct, cone, cyl, point);
  return (color);
}

int                     main(int argc, char **argv)
{
  t_main *main;

  main = (t_main *)malloc(sizeof(t_main));
  main->img = (t_img *)malloc(sizeof(t_img) * 1);
  main->mlx = mlx_init();
  main->win = mlx_new_window(main->mlx, WIDTH, HEIGHT, "RayTracer V1");
  main->img->img_ptr = mlx_new_image(main->mlx, WIDTH, HEIGHT);
  parcour_all(main);
  mlx_put_image_to_window(main->mlx, main->win, main->img->img_ptr, 0, 0);
  mlx_key_hook(main->win, event_mlx, main);
  mlx_loop(main->mlx);
  return (0);
}
