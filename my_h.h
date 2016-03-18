#ifndef MY_H
#define MY_H
# include <mlx.h>
# include <stdio.h>
# include <stdlib.h>
# include <unistd.h>
# include <math.h>
# define RGB(r, g, b)(256 * 256 * (int)(r) + 256 * (int)(g) + (int)(b))
# define WIDTH 800
# define HEIGHT 600

typedef struct                                                                  s_img
{
  double cos;
  void                            *mlx;
  void                            *win;
  unsigned long           img_color;
  char                            *data;
  void                            *img_ptr;
  int                                     sizeline;
  int                                     bpp;
  int                                     x;
  int                                     y;
  int                                     endian;
  void                            *mlx_ptr;
}                                                                                               t_img;

typedef struct  s_pos
{
  float         k_plan;
  float         k_sph;
  float         k_cyl;
  float         k_cone;
  float         xp;
  float         yp;
  float         zp;
}               t_pos;



typedef struct  s_coord
{
  float         vx;
  float         vy;
  float         vz;
  float         rayon;
  float         xeye;
  float         yeye;
  float         zeye;
  float         angle;
}               t_coord;

typedef struct  s_lum
{
  float         x;
  float         y;
  float         z;
  float         vx;
  float         vy;
  float         vz;
  float         nx;
  float         ny;
  float         nz;
}               t_lum;

typedef struct  s_struct_point
{
  t_coord       *coord;
  t_coord       *sph;
  t_coord       *cyl;
  t_coord       *cone;
  t_pos         *pos;
  t_lum         *lum;
  t_img 	*img;
  int count;
}               t_point;



typedef struct s_vec3d
{
  double x;
  double y;
  double z;
  double v;
  double rayon;
  double         xp;
  double         yp;
  double        zp;
}		t_vec3d;


typedef struct s_sphere
{
  t_vec3d *pos;
  double xt;
  double yt;
  double zt;
  double radius;
}		t_sphere;


typedef struct s_ray
{
  t_vec3d *o;
  t_vec3d *d;
  double vx;
  double vy;
  double vz;
  double xeye;
  double yeye;
  double zeye;
}		t_ray;


typedef struct s_light
{
  double x;
  double y;
  double z;
  double dx;
  double dy;
  double dz;
  double ox;
  double oy;
  double oz;
}		t_light;

typedef struct                                                                  s_struct
{
  int                                     fracnum;
  void                            *mlx;
  void                            *win;
  t_img                           *img;
  int color;
  double decx;
  double decy;
  double decz;
  t_sphere *spheres;
  t_sphere *cones;
  t_sphere *cyls;
  t_vec3d *plan;
  double xp;
  double yp;
  double zp;
  t_ray *rays;
  t_light *lights;
  t_point *point;
  int count;
}                                                                                               t_struct;

typedef struct s_main
{
  void                            *mlx;
  void                            *win;
  t_img                           *img;
  t_vec3d 	*sphere;
  t_vec3d 	*cyl;
  t_vec3d 	*cone;
  t_vec3d *plan;
  t_coord       *eye;
  t_lum         *lum;
  int color;
}	t_main;
int definee_color(int x, int y, t_main *main);
void give_value(int x, int y, t_sphere *sphere, t_ray *ray, t_light *light, t_struct *mystruct, t_sphere *cone, t_sphere *cyl);
int             my_pixel_put_to_image2(t_img *img, int x, int y, int color);
void give_v(t_vec3d *plan, t_sphere *sphere, t_sphere *cone, t_sphere *cyl, double *v_value, t_point *point);
int light_plan(t_point *point, t_struct *mystruct);
int plan_int(t_struct *mystruct);
void sphere_int(t_sphere *sphere, t_ray *ray);
void cone_int(t_sphere *cone, t_ray *ray);
void cyl_int(t_sphere *cyl, t_ray *ray);
void parcour_all(t_main *main);
int define_color(int x, int y, t_struct *mystruct);
void give_cylinder(t_sphere *cyl, t_light *light, t_ray *ray, double v);
void give_cone(t_sphere *cone, t_light *light, t_ray *ray, double v);
void give_sphere(t_sphere *sphere, t_light *light, t_ray *ray, double v, t_struct *mystruct);
int                     event_mlx(int keycode, t_main *main);
int shadow_cone(t_sphere *cone, t_light *light, double v, t_ray *ray);
int shadow_cyl(t_sphere *cyl, t_light *light, double v, t_ray *ray);
int shadow_sphere(t_sphere *sphere, t_light *light, t_ray *ray, t_struct *mystruct, double v);
double mypow(int a, int n);
void my_pixel_put_to_image(t_img *myimg, int x, int y, int color);
#endif
