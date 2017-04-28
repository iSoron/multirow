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

#define _XOPEN_SOURCE 700

#include <math.h>
#include <stdio.h>

#include <multirow/double.h>
#include <multirow/util.h>
#include <multirow/lp.h>
#include <multirow/lfree2d.h>

int LFREE_2D_init(struct LFreeSet2D *set,
                  int n_vertices,
                  int n_lattice_points,
                  int n_halfspaces)
{
    int rval = 0;

    set->vertices = (double *) malloc(2 * n_vertices * sizeof(double));
    set->lattice_points = (double *) malloc( 2 * n_lattice_points * sizeof(double));
    set->halfspaces = (double *) malloc(2 * n_halfspaces * sizeof(double));

    abort_if(!set->vertices, "could not allocate set->vertices");
    abort_if(!set->lattice_points, "could not allocate set->lattice_points");
    abort_if(!set->halfspaces, "could not allocate set->halfspaces");

CLEANUP:
    return rval;
}

void LFREE_2D_free(struct LFreeSet2D *set)
{
    if(set->vertices) free(set->vertices);
    if(set->lattice_points) free(set->lattice_points);
    if(set->halfspaces) free(set->halfspaces);
}

/**
 * Computes the bounding box around the given set.
 *
 * @param[in]  set  the set whose bounding box should be computed
 * @param[out] lb   two-dimensional vector representing lower-left corner
 *                  of the bounding box
 * @param[out] ub   two-dimensional vector representing upper-right corner
 *                  of the bounding box
 */
int LFREE_2D_get_bounding_box(const struct LFreeSet2D *set,
                              int *lb,
                              int *ub)
{
    int rval = 0;

    ub[0] = ub[1] = INT_MIN;
    lb[0] = lb[1] = INT_MAX;

    for (int i = 0; i < set->n_vertices; i++)
    {
        double *vertex = &set->vertices[2 * i];

        ub[0] = (int) fmax(ub[0], ceil(vertex[0]));
        ub[1] = (int) fmax(ub[1], ceil(vertex[1]));

        lb[0] = (int) fmin(lb[0], floor(vertex[0]));
        lb[1] = (int) fmin(lb[1], floor(vertex[1]));
    }

CLEANUP:
    return rval;
}

int LFREE_2D_print_set(const struct LFreeSet2D *set)
{
    int rval = 0;

    log_debug("  f=%12.6lf %12.6lf\n", set->f[0], set->f[1]);

    for (int i = 0; i < set->n_vertices; i++)
    {
        double *vertex = &set->vertices[2 * i];
        log_debug("  v%-3d=%12.6lf %12.6lf\n", i, vertex[0], vertex[1]);
    }

    for (int i = 0; i < set->n_lattice_points; i++)
    {
        double *lattice_point = &set->lattice_points[2 * i];
        log_debug("  t%-3d=%12.6lf %12.6lf\n", i, lattice_point[0], lattice_point[1]);
    }

    for(int i = 0; i < set->n_halfspaces; i++)
    {
        double *halfspace = &set->halfspaces[2 * i];
        log_debug("  h%-3d=%12.6lf %12.6lf\n", i, halfspace[0], halfspace[1]);
    }

CLEANUP:
    return rval;
}

/**
 * Computes the halfspace representation for a given set.
 *
 * To construct a LFreeSet2D, it is only necessary to specify
 * the vertices and the interior point. This function computes
 * the halfspace representation. It assumes that the vertices are
 * sorted in either clockwise or counter-clockwise order.
 * The set is modified in-place.
 *
 * @param set  the set whose halfspace representation should
 *             be computed
 *
 */
int LFREE_2D_compute_halfspaces(struct LFreeSet2D *set)
{
    int rval = 0;

    int k = 0;
    int expected_sgn = 0;
    double *f = set->f;

    for (int i0 = 0; i0 < set->n_vertices; i0++)
    {
        int i1 = (i0 + 1) % set->n_vertices;
        int i2 = (i0 + 2) % set->n_vertices;

        double *v0 = &set->vertices[2 * i0];
        double *v1 = &set->vertices[2 * i1];
        double *v2 = &set->vertices[2 * i2];

        double x0 = v0[0], y0 = v0[1];
        double x1 = v1[0], y1 = v1[1];
        double x2 = v2[0], y2 = v2[1];

        int sgn = DOUBLE_sgn(
                -x1 * y0 + x2 * y0 + x0 * y1 - x2 * y1 - x0 * y2 + x1 * y2);

        abort_if(sgn == 0, "vertices should not be aligned");
        if(expected_sgn == 0) expected_sgn = sgn;
        abort_if(expected_sgn != sgn, "wrong orientation of vertices");

        double *halfspace = &set->halfspaces[2 * (k++)];

        halfspace[0] = y0 - y1;
        halfspace[1] = x1 - x0;

        double rhs = y0 * x1 - x0 * y1;
        rhs -= (y0 - y1) * f[0];
        rhs -= (x1 - x0) * f[1];

        halfspace[0] /= rhs;
        halfspace[1] /= rhs;
    }

    set->n_halfspaces = k;

CLEANUP:
    return rval;
}

/**
 * Translates a given set by a certain amount.
 *
 * @param set  the set to be translated
 * @param dx   offset to be applied to the first component
 * @param dy   offset to be applied to the second component
 */
int LFREE_2D_translate_set(struct LFreeSet2D *set, double dx, double dy)
{
    set->f[0] += dx;
    set->f[1] += dy;

    for (int i = 0; i < set->n_vertices; i++)
    {
        double *vertex = &set->vertices[2 * i];
        vertex[0] += dx;
        vertex[1] += dy;
    }

    for (int i = 0; i < set->n_lattice_points; i++)
    {
        double *lattice_point = &set->lattice_points[2 * i];
        lattice_point[0] += dx;
        lattice_point[1] += dy;
    }

    return 0;
}

/**
 * Computes the inverse of the given matrix and left-multiplies it
 * by the given ray. The matrix is given by [[x1,y1],[x2,y2]]. The
 * ray is modified in-place.
 */
static void apply_inverse(double x1,
                          double y1,
                          double x2,
                          double y2,
                          double *ray)
{
    double a = ray[0];
    double b = ray[1];
    double det = x2 * y1 - x1 * y2;

    ray[0] = (b * x2 - a * y2) / det;
    ray[1] = (a * y1 - b * x1) / det;
}

/**
 * Applies a unimodular affine transformation to the set, making
 * lattice points 1 and 2 orthogonal. Returns the transformation
 * applied.
 *
 * @param[in] set the set to be transformed
 * @param[out] m two-by-two matrix representing the transformation
 *               applied
 */
static int make_t1_t2_orthogonal(struct LFreeSet2D *set, double *m)
{
    int rval = 0;

    double x1 = set->lattice_points[2];
    double y1 = set->lattice_points[3];

    double x2 = set->lattice_points[4];
    double y2 = set->lattice_points[5];

    double det = x2 * y1 - x1 * y2;
    abort_if(DOUBLE_iszero(det), "determinant should not be zero");

    apply_inverse(x1, y1, x2, y2, set->f);

    for (int i = 0; i < set->n_vertices; i++)
    {
        double *vertex = &set->vertices[2 * i];
        apply_inverse(x1, y1, x2, y2, vertex);
    }

    for (int i = 0; i < set->n_lattice_points; i++)
    {
        double *lattice_point = &set->lattice_points[2 * i];
        apply_inverse(x1, y1, x2, y2, lattice_point);
    }

    apply_inverse(x1, y1, x2, y2, &m[0]);
    apply_inverse(x1, y1, x2, y2, &m[2]);

CLEANUP:
    return rval;
}

static int apply_m_to_vect(double *v, const double *m)
{
    double v0 = v[0], v1 = v[1];
    v[0] = m[0] * v0 + m[1] * v1;
    v[1] = m[2] * v0 + m[3] * v1;

    return 0;
}

/**
 * Applies a unimodular affine transformation to the set, making
 * the four lattice points the unit square.
 *
 * Must be called on quadrilaterals only, and assumes that the three
 * first lattice points already correspond to (0,0), (0,1) and (1,0).
 * The set is modified in place. The transformation is also applied
 * to a given matrix pre_m.
 */
static int make_square(struct LFreeSet2D *set, double *pre_m)
{
    int rval = 0;
    double m[4];

    abort_if(set->n_vertices != 4, "make_square requires quadrilaterals");

    double *t = &set->lattice_points[2 * 3];
    if(t[0] == 1 && t[1] == 1)
    {
        m[0] = 1; m[1] = 0;
        m[2] = 0; m[3] = 1;
    }
    else if(t[0] == 1 && t[1] == -1)
    {
        m[0] = 1; m[1] = 0;
        m[2] = 1; m[3] = 1;
    }
    else if(t[0] == -1 && t[1] == 1)
    {
        m[0] = 1; m[1] = 1;
        m[2] = 0; m[3] = 1;
    }
    else
    {
        abort_if(1, "invalid state");
    }

    rval = LFREE_2D_transform_set(set, m);
    abort_if(rval, "LFREE_2D_transform_set failed");

    apply_m_to_vect(&pre_m[0], m);
    apply_m_to_vect(&pre_m[2], m);

CLEANUP:
    return rval;
}


int find_thin_direction(const struct LFreeSet2D *set,
                        double *thin_direction,
                        double *width)
{
    int rval = 0;

    int n_candidates = 3;
    double candidates[] = {
            0,  1,
            1,  0,
            1,  1
    };

    double best_width = INFINITY;
    double *best_candidate = 0;

    for (int i = 0; i < n_candidates; i++)
    {
        double *d = &candidates[2 * i];
        double max_val = -INFINITY;
        double min_val = INFINITY;

        for (int j = 0; j < set->n_vertices; j++)
        {
            double *vertex = &set->vertices[2 * j];
            double val = vertex[0] * d[0] + vertex[1] * d[1];

            max_val = fmax(max_val, val);
            min_val = fmin(min_val, val);
        }

        double delta = fabs(max_val - min_val);

        log_verbose("d=%.3lf %.3lf width=%.3lf\n", d[0], d[1], delta);

        if(delta < best_width)
        {
            best_width = delta;
            best_candidate = d;
        }
    }

    memcpy(thin_direction, best_candidate, 2 * sizeof(double));
    *width = best_width;

CLEANUP:
    return rval;
}

int LFREE_2D_transform_set(struct LFreeSet2D *set, const double *m)
{
    int rval = 0;

    rval = apply_m_to_vect(set->f, m);
    abort_if(rval, "apply_m_to_vect failed");

    for (int i = 0; i < set->n_vertices; i++)
    {
        double *vertex = &set->vertices[2 * i];

        rval = apply_m_to_vect(vertex, m);
        abort_if(rval, "apply_m_to_vect failed");
    }

    for (int i = 0; i < set->n_lattice_points; i++)
    {
        double *lattice_point = &set->lattice_points[2 * i];

        rval = apply_m_to_vect(lattice_point, m);
        abort_if(rval, "apply_m_to_vect failed");
    }

CLEANUP:
    return rval;
}

int LFREE_2D_preprocess(struct LFreeSet2D *set, double *pre_m, double *center)
{
    int rval = 0;
    double direction[2], m[4], width;
    double translate[2];

    pre_m[0] = 1; pre_m[1] = 0;
    pre_m[2] = 0; pre_m[3] = 1;

    rval = LFREE_2D_translate_set(set, -set->lattice_points[0],
            -set->lattice_points[1]);
    abort_if(rval, "LFREE_2D_translate set failed");

    rval = make_t1_t2_orthogonal(set, pre_m);
    abort_if(rval, "make_t1_t2_orthogonal failed");

    if(set->n_vertices == 4)
    {
        rval = make_square(set, pre_m);
        abort_if(rval, "make_square failed");
    }

    rval = find_thin_direction(set, direction, &width);
    abort_if(rval, "find_thin_direction failed");

    log_debug("d=%.2lf %.2lf\n", direction[0], direction[1]);

    if(direction[0] == 0 && direction[1] == 1)
    {
        m[0] = 1; m[1] = 0;
        m[2] = 0; m[3] = 1;
        translate[0] = translate[1] = 0;
    }
    else if(direction[0] == 1 && direction[1] == 0)
    {
        m[0] = 0; m[1] = 1;
        m[2] = 1; m[3] = 0;
        translate[0] = translate[1] = 0;
    }
    else if(direction[0] == 1 && direction[1] == 1)
    {
        m[0] =  1; m[1] =  0;
        m[2] =  -1; m[3] =  -1;
        translate[0] = 0;
        translate[1] = 1;
    }
    else
    {
        abort_if(1, "invalid direction");
    }

    apply_m_to_vect(&pre_m[0], m);
    apply_m_to_vect(&pre_m[2], m);

    rval = LFREE_2D_transform_set(set, m);
    abort_if(rval, "LFREE_2D_transform_set failed");

    rval = LFREE_2D_translate_set(set, translate[0], translate[1]);
    abort_if(rval, "LFREE_2D_translate_set failed");

    double max_v1 = -INFINITY;
    double min_v1 = INFINITY;

    for(int i = 0; i < set->n_vertices; i++)
    {
        double *v = &set->vertices[2 * i];
        max_v1 = fmax(max_v1, v[1]);
        min_v1 = fmin(min_v1, v[1]);
    }

    *center = (max_v1 + min_v1) / 2;
    width = max_v1 - min_v1;

    log_debug("width = %.2lf\n", width);
    abort_iff(width >= 3, "lattice width too large: %.2lf", width);

CLEANUP:
    return rval;
}

int LFREE_2D_read_next(FILE *file, struct LFreeSet2D *set)
{
    int rval = 0;
    ssize_t count;

    double *f = set->f;
    count = fscanf(file, "%lf %lf ", &f[0], &f[1]);
    abort_if(count != 2, "could not read f");

    count = fscanf(file, "%d ", &set->n_vertices);
    abort_if(count != 1, "could not read n_vertices");
    abort_if(set->n_vertices != 3 && set->n_vertices != 4,
            "only triangles or quadrilaterals supported");

    for(int i = 0; i < set->n_vertices; i++)
    {
        double *vertex = &set->vertices[2 * i];
        count = fscanf(file, "%lf %lf ", &vertex[0], &vertex[1]);
        abort_iff(count != 2, "could not read vertex %d", i+1);
    }

    count = fscanf(file, "%d ", &set->n_lattice_points);
    abort_if(count != 1, "could not read n_lattice_points");
    abort_if(set->n_lattice_points != set->n_vertices, "number of vertices "
            "should be the same as number of lattice points");

    for(int i = 0; i < set->n_lattice_points; i++)
    {
        double *lattice_point = &set->lattice_points[2 * i];
        count = fscanf(file, "%lf %lf ", &lattice_point[0], &lattice_point[1]);
        abort_iff(count != 2, "could not read lattice_point %d", i+1);
    }

CLEANUP:
    return rval;
}

