#include "my_h.h"

int                     event_mlx(int keycode, t_main *main)
{

  if (keycode == 53)
    exit(1);
  if(keycode == 123)
    {
      //gauche
      /*      mlx_destroy_image(mystruct->mlx, mystruct->img->img_ptr);
      mystruct->img->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
      mystruct->count++;
      parcour_all(mystruct);
      mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);*/
    }
  if(keycode == 124)
    {
      //droite
    }
  return (0);
}

double mypow(int a, int n)
{
  double res;
  int i;

  i = 1;
  res = a;
  while(i < n)
    {
      res = res * a;
      i++;
    }
  return (res);
}
