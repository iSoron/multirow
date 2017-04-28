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

#ifndef LIFTING_H
#define LIFTING_H

#include <multirow/lfree2d.h>

int LIFTING_2D_psi(int n_halfspaces,
                   const double *halfspaces,
                   const double *ray,
                   double *value);

int LIFTING_2D_optimize_continuous(int n_halfspaces,
                                   const double *halfspaces,
                                   double alpha2,
                                   double *alpha1,
                                   double *value);

int LIFTING_2D_lift_fixed(int n_halfspaces,
                          const double *halfspaces,
                          const double *ray,
                          double k1,
                          double *opt);

int LIFTING_2D_naive(int n_halfspaces,
                     const double *halfspaces,
                     const double *ray,
                     const int *lb,
                     const int *ub,
                     double *value);

int LIFTING_2D_bound(int n_halfspaces,
                     const double *halfspaces,
                     const double *ray,
                     double *value);

int LIFTING_2D_verify(struct LFreeSet2D *set);

#endif //LIFTING_H
