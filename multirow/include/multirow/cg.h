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

#ifndef MULTIROW_CG_H
#define MULTIROW_CG_H

#include "lp.h"

struct CG
{
    struct LP *lp;
    struct Row **tableau_rows;
    int *cstat;
    int *rstat;
    double *ub;
    double *lb;
    int nrows;
    int ncols;
    char *column_types;

    double last_obj_value;

    double *integral_solution;
    double *basic_solution;
    double *current_solution;
};

typedef int (*SingleRowGeneratorCallback)(const struct Row *row,
                                          char *column_types,
                                          struct Row *cut);

typedef int (*MultirowGeneratorCallback)(int nrows,
                                         struct Row **rows,
                                         char *column_types,
                                         struct Row *cut);

int CG_init(struct LP *lp,
            char *column_types,
            struct CG *cg);

void CG_free(struct CG *cg);

int CG_add_single_row_cuts(struct CG *cg,
                           SingleRowGeneratorCallback generate);

int CG_add_multirow_cuts(struct CG *cg,
                         int nrows,
                         MultirowGeneratorCallback generate);

int CG_set_integral_solution(struct CG *cg,
                             double *valid_solution);

int CG_set_basic_solution(struct CG *cg,
                          double *basic_solution);

int CG_extract_rays_from_rows(int nrows,
                              struct Row **rows,
                              double *rays,
                              int *nrays,
                              int *variable_to_ray,
                              double *ray_scale,
                              int *indices,
                              int *nz);

int CG_boost_variable(int var,
                      double factor,
                      int nrows,
                      double *rays,
                      int *variable_to_ray,
                      double *ray_scale,
                      int *indices,
                      int nz);

int CG_find_ray(int dim,
                const double *rays,
                int nrays,
                const double *r,
                int *found,
                double *scale,
                int *index);

#endif //MULTIROW_CG_H
