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

#ifndef MULTIROW_GREEDY_2D_H
#define MULTIROW_GREEDY_2D_H

int GREEDY_2D_bound(const double *rays,
                    const double *bounds,
                    int nrays,
                    const double *f,
                    const double *p,
                    double *epsilon,
                    double *v1,
                    double *v2,
                    int *index1,
                    int *index2);

int GREEDY_2D_generate_cut(const double *rays,
                           int nrays,
                           const double *f,
                           double *bounds);


#endif //MULTIROW_GREEDY_2D_H
