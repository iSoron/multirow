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

#ifndef MULTIROW_DOUBLE_H
#define MULTIROW_DOUBLE_H

#include <multirow/rational.h>

int DOUBLE_sgn(double a);

int DOUBLE_cmp(double a,
               double b);

int DOUBLE_geq(double a,
               double b);

int DOUBLE_leq(double a,
               double b);

int DOUBLE_eq(double a,
              double b);

int DOUBLE_neq(double a,
               double b);

double DOUBLE_max(double a,
                  double b);

int DOUBLE_to_rational(double a,
                       long max_den,
                       Rational r);

#define DOUBLE_iszero(a) (fabs(a) < EPSILON)

double DOUBLE_random(double min, double max);

int DOUBLE_sum(const double *list, const int count, double *result);

#endif //MULTIROW_DOUBLE_H
