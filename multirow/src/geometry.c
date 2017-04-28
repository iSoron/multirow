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

#define _GNU_SOURCE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <multirow/geometry.h>
#include <multirow/double.h>
#include <multirow/util.h>

struct SortPair
{
    int index;
    void *data;
};

#if defined(__APPLE__)
static int _qsort_cmp_angle_distance(void *p_orig,
                                     const void *p1,
                                     const void *p2)
#else
static int _qsort_cmp_angle_distance(const void *p1,
                                     const void *p2,
                                     void *p_orig)
#endif
{
    double *orig = (double *) p_orig;
    double *r1 = (double*) p1;
    double *r2 = (double*) p2;

    double norm1 = fabs(r1[0] - orig[0]) + fabs(r1[1] - orig[1]);
    double norm2 = fabs(r2[0] - orig[0]) + fabs(r2[1] - orig[1]);

    double angle1 = atan2(r1[0] - orig[0], r1[1] - orig[1]);
    double angle2 = atan2(r2[0] - orig[0], r2[1] - orig[1]);

    if(DOUBLE_neq(angle1, angle2))
        return sign(angle2 - angle1);
    else
        return sign(norm1 - norm2);
}

static double turn_direction(const double *r1,
                      const double *r2,
                      const double *r3)
{
    return (r2[0]-r1[0]) * (r3[1]-r1[1]) - (r2[1]-r1[1]) * (r3[0]-r1[0]);
}

int chull_2d(const double *original_points,
             int npoints,
             double *vertices,
             int *nvertices)
{
    int rval = 0;
    int pivot_idx;
    double pivot[2];

    int *stack = 0;
    int stack_length = 0;
    double *points = 0;

    stack = (int*) malloc(npoints * sizeof(int));
    points = (double*) malloc(2 * npoints * sizeof(double));
    abort_if(!stack, "could not allocate stack");
    abort_if(!points, "could not allocate points");

    memcpy(points, original_points, 2 * npoints * sizeof(double));

    // find point with lowest y-coordinate
    pivot_idx = 0;
    memcpy(pivot, points, 2 * sizeof(double));

    for (int i = 0; i < npoints; i++)
    {
        double *point = &points[2*i];

        if(point[1] > pivot[1]) continue;
        if(point[1] == pivot[1] && point[0] > pivot[0]) continue;

        memcpy(pivot, point, 2 * sizeof(double));
        pivot_idx = i;
    }

    // remove pivot
    swap(points[2*(npoints-1)], points[2*pivot_idx], double);
    swap(points[2*(npoints-1)+1], points[2*pivot_idx+1], double);
    npoints--;

    // sort points counterclockwise
#if defined(__APPLE__)
    qsort_r(points, npoints, 2 * sizeof(double), pivot, _qsort_cmp_angle_distance);
#else
    qsort_r(points, npoints, 2 * sizeof(double), _qsort_cmp_angle_distance,
            pivot);
#endif
    
    // adds pivot and first point to the stack
    stack[stack_length++] = npoints;
    stack[stack_length++] = 0;

    int k = 1;
    while(k < npoints)
    {
        double *p0 = &points[2*stack[stack_length-2]];
        double *p1 = &points[2*stack[stack_length-1]];
        double *p2 = &points[2*k];

        double t = turn_direction(p0, p1, p2);

        if(DOUBLE_leq(t, 0))
            stack_length--;

        if(DOUBLE_geq(t, 0))
        {
            stack[stack_length++] = k;
            k++;
        }
    }

    for(int i = 0; i < stack_length; i++)
    {
        int j = stack_length - i - 1;
        vertices[2*j] = points[2*stack[i]];
        vertices[2*j+1] = points[2*stack[i]+1];
    }

    *nvertices = stack_length;

    CLEANUP:
    if(points) free(points);
    if(stack) free(stack);
    return rval;
}

/*
 * Find lambda such that lambda * p lies on the line connecting ray1 to ray2
 */
int scale_vector_to_line(const double *ray1,
                         const double *ray2,
                         const double *p,
                         double *lambda)
{
    int rval = 0;

    double a = ray1[0], b = ray1[1];
    double c = ray2[0], d = ray2[1];

    double norm1 = fabs(ray1[0]) + fabs(ray1[1]);
    double norm2 = fabs(ray2[0]) + fabs(ray2[1]);

    log_verbose("r1=%12.8e %12.8e\n", ray1[0], ray1[1]);
    log_verbose("r2=%12.8e %12.8e\n", ray2[0], ray2[1]);
    log_verbose("p=%12.8e %12.8e\n\n", p[0], p[1]);

    if (norm1 > norm2)
    {
        swap(a, c, double);
        swap(b, d, double);
        swap(norm1, norm2, double);
    }

    if (DOUBLE_iszero(norm1 / norm2))
    {
        *lambda = (b * c - a * d) / (c * p[1] - d * p[0]);
    }
    else
    {
        *lambda = b * c - a * d;
        *lambda /= p[0] * (b - d) - p[1] * (a - c);
    }

CLEANUP:
    return rval;
}
