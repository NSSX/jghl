#include "my_h.h"

void parcour_all(t_main *main)
{
  int i;
  int j;
  int mypixelcolor;

  mypixelcolor = 0;
  j = 0;
  while(j < HEIGHT)
    {
      i = 0;
      while(i < WIDTH)
	{
          mypixelcolor = definee_color(i, j, main);
          if(mypixelcolor != 0)
            {
	      my_pixel_put_to_image2(main->img, i, j, mypixelcolor);
            }
          i++;
	}
      j++;
    }
}

