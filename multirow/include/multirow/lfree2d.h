/* Copyright (c) 2015 Alinson Xavier
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LFREE_2D_H
#define LFREE_2D_H

#include <stdio.h>

struct LFreeSet2D
{
    double f[2];

    int n_vertices;
    double *vertices;

    int n_lattice_points;
    double *lattice_points;

    int n_halfspaces;
    double *halfspaces;
};

struct RayList
{
    double *values;
    int nrays;
    int dim;
};

struct ConvLFreeSet
{
    double *f;
    struct RayList rays;
    int nrows;
    double *beta;
};

int LFREE_2D_init(struct LFreeSet2D *set,
                  int n_vertices,
                  int n_lattice_points,
                  int max_n_halfspaces);

void LFREE_2D_free(struct LFreeSet2D *set);

int LFREE_2D_read_next(FILE *file, struct LFreeSet2D *set);

int LFREE_2D_compute_halfspaces(struct LFreeSet2D *set);

int LFREE_2D_preprocess(struct LFreeSet2D *set, double *m, double *center);

int LFREE_2D_transform_set(struct LFreeSet2D *set, const double *m);

int LFREE_2D_print_set(const struct LFreeSet2D *set);

int LFREE_2D_translate_set(struct LFreeSet2D *set, double dx, double dy);

int LFREE_2D_get_bounding_box(const struct LFreeSet2D *set, int *lb, int *ub);

void LFREE_push_ray(struct RayList *list, const double *ray);

double* LFREE_get_ray(const struct RayList *list, int index);

void LFREE_free_ray_list(struct RayList *list);

int LFREE_init_ray_list(struct RayList *list, int dim, int capacity);

int LFREE_init_conv(struct ConvLFreeSet *lfree, int dim, int max_nrays);

void LFREE_free_conv(struct ConvLFreeSet *lfree);

int LFREE_print_set(const struct ConvLFreeSet *lfree);

#endif //LFREE_2D_H
