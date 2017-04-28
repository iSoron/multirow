/* Copyright (c) 2016 Laurent Poirrier
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

#ifndef LIFTING_MIP_H
#define LIFTING_MIP_H

extern double MIP_TIME_OPTIMIZE;
extern double MIP_TIME_CREATE;

int LIFTING_2D_mip_init();

void LIFTING_2D_mip_cleanup();

int LIFTING_2D_mip(int n_halfspaces,
                   const double *halfspaces,
                   const double *ray,
                   double *value);

#endif
