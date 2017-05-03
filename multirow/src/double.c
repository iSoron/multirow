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

#include <math.h>
#include <float.h>
#include <stdlib.h>

#include <multirow/double.h>
#include <multirow/util.h>

int DOUBLE_sgn(double a)
{
    return (a > 0 ? 1 : (a < 0 ? -1 : 0));
}

int DOUBLE_cmp(double a,
               double b)
{
    if(isnan(a) || isnan(b)) return 1;

    double size = min(fabs(a) + fabs(b), DBL_MAX);
    double zero_cutoff = EPSILON * size;
    double diff = (a-b);

    // Do not use relative epsilon for numbers close to zero
    if(fabs(a) < EPSILON || fabs(b) < EPSILON) zero_cutoff = EPSILON;

    return (fabs(diff) < zero_cutoff ? 0 : DOUBLE_sgn(diff));
}

int DOUBLE_geq(double a,
               double b)
{
    return DOUBLE_cmp(a, b) >= 0;
}

int DOUBLE_leq(double a,
               double b)
{
    return DOUBLE_cmp(a, b) <= 0;
}

int DOUBLE_eq(double a,
              double b)
{
    return DOUBLE_cmp(a, b) == 0;
}

int DOUBLE_neq(double a,
               double b)
{
    return DOUBLE_cmp(a, b) != 0;
}

double DOUBLE_max(double a,
                  double b)
{
    return (DOUBLE_cmp(a,b) > 0 ? a : b);
}

double DOUBLE_random(double min, double max)
{
    return min + ((double) rand() / ((double) RAND_MAX / (max - min)));
}

int DOUBLE_to_rational(double a,
                       long max_den,
                       Rational r)
{
    int rval = 0;

    abort_if(!isfinite(a), "a must be finite");
    abort_if(max_den < 1, "max_den must be positive");

    long shift = floor(a);
    a = frac(a);

    long n1 = 0, d1 = 1;
    long n2 = 1, d2 = 1;

    log_verbose("a=%.8lf\n", a);

    while(d1 + d2 < max_den)
    {
        double mid = ((double) (n1 + n2) / (d1 + d2));
        log_verbose("  d1=%ld %d\n", n1, d1);
        log_verbose("  d2=%ld %d\n", n2, d2);
        log_verbose(" mid=%.8lf\n", mid);


        if(a > mid)
        {
            n1 = n1 + n2;
            d1 = d1 + d2;
        }
        else
        {
            n2 = n1 + n2;
            d2 = d1 + d2;
        }
    }

    double err1 = fabs(a - ((double) n1 / d1));
    double err2 = fabs(a - ((double) n2 / d2));

    if(err1 < err2)
    {
        r->num = n1;
        r->den = d1;
    }
    else
    {
        r->num = n2;
        r->den = d2;
    }

    r->num += ((r->den) * shift);

CLEANUP:
    return rval;
}

static int _qsort_double_cmp(const void *a, const void *b)
{
    double ia = *((double*) a);
    double ib = *((double*) b);
    if (ia > ib) return -1;
    if (ia < ib) return 1;
    return 0;
}

int DOUBLE_sum(const double *list, const int count, double *result)
{
    int rval = 0;

    double *list_copy = 0;
    list_copy = (double*) malloc(count * sizeof(double));
    abort_if(!list_copy, "could not allocate list_copy");

    memcpy(list_copy, list, count * sizeof(double));
    qsort(list_copy, (size_t) count, sizeof(double), _qsort_double_cmp);

    *result = 0;
    for(int i = 0; i < count; i++)
        *result += list_copy[i];

CLEANUP:
    if(list_copy) free(list_copy);
    return rval;
}
