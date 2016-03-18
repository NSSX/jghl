/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   all_give.c                                         :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: avella <marvin@42.fr>                      +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2016/03/17 15:48:12 by avella            #+#    #+#             */
/*   Updated: 2016/03/17 15:48:33 by avella           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "my_h.h"

void	give_cylinder(t_sphere *cyl, t_light *light, t_ray *ray, double v)
{
	cyl->xt =ray->o->x + (cyl->pos->x *v);
	cyl->yt= ray->o->y + (cyl->pos->y * v);
	light->dx = light->x - cyl->xt;
	light->dy = light->y - cyl->yt;
}


void	give_cone(t_sphere *cone, t_light *light, t_ray *ray, double v)
{
	cone->xt = ray->o->x + (cone->pos->x * v);
	cone->yt = ray->o->y + (cone->pos->y * v);
	cone->zt = ray->o->z + (cone->pos->z * v);
	light->dx = light->x - cone->xt;
	light->dy = light->y - cone->yt;
	light->dz = light->z - cone->zt;
}


void	give_sphere(t_sphere *sphere, t_light *light, t_ray *ray, double v, t_struct *mystruct)
{
	mystruct->xp = ray->o->x + (sphere->pos->x * v);
	mystruct->yp = ray->o->y + (sphere->pos->y * v);
	mystruct->zp = ray->o->z + (sphere->pos->z * v);
	mystruct->lights->dx = mystruct->lights->x - mystruct->xp;
	mystruct->lights->dy = mystruct->lights->y - mystruct->yp;
	mystruct->lights->dz = mystruct->lights->z - mystruct->zp;
}


void	give_value(int x, int y, t_sphere *sphere, t_ray *ray, t_light *light, t_struct *mystruct, t_sphere *cone, t_sphere *cyl)
{

	mystruct->lights->ox = 0;
	mystruct->lights->oy = 0;
	mystruct->lights->oz = 100;
	light->x = -3000;
	light->y = 0;
	light->z = 2000;
	ray->o->x = -3000;
	ray->o->y = 0;
	ray->o->z = 500;
	ray->d->x = 0;
	ray->d->y = (WIDTH / 2) - x;
	ray->d->z = (HEIGHT / 2) - y;
	sphere->pos->x = ray->d->x + mystruct->decx - 400;
	sphere->pos->y = ray->d->y + mystruct->decy;
	sphere->pos->z = ray->d->z + mystruct->decz;
	sphere->radius = 800;
	cone->pos->x = ray->d->x + mystruct->decx;
	cone->pos->y = ray->d->y + mystruct->decy + 100;
	cone->pos->z = ray->d->z + mystruct->decz;
	cyl->pos->x = ray->d->x + mystruct->decx;
	cyl->pos->y = ray->d->y + mystruct->decy - 100;
	cyl->pos->z = ray->d->z + mystruct->decz;
	cyl->radius = 150;

}

void	give_v(t_vec3d *plan, t_sphere *sphere, t_sphere *cone, t_sphere *cyl, double *v_value, t_point *point)
{
	v_value[0] = plan->v;//point->pos->k_plan;//plan->v;
	v_value[1] = cone->pos->v;
	v_value[2] = cyl->pos->v;
	v_value[3] = sphere->pos->v;
}
