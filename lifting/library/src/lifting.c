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
#include <lifting/lifting.h>

/**
 * Verifies if the set is well formed.
 */
int LIFTING_2D_verify(struct LFreeSet2D *set)
{
    int rval = 0;

    abort_if(set->n_halfspaces == 0, "Halfspaces not found");

    for(int i = 0; i < set->n_lattice_points; i++)
    {
        double *t = &set->lattice_points[2 * i];
        double r[2] = {t[0] - set->f[0],
                       t[1] - set->f[1]};

        double value;

        rval = LIFTING_2D_psi(set->n_halfspaces, set->halfspaces, r, &value);
        abort_if(rval, "LIFTING_2D_psi failed");

        double delta = fabs(value - 1);
        if(delta > 0.0001)
        {
            log_debug("Lattice point (%.2lf, %.2lf) is "
                "not on the boundary (delta=%.6lf)", t[0], t[1], delta);
            rval = 1;
            goto CLEANUP;
        }
    }

CLEANUP:
    return rval;
}

int LIFTING_2D_psi(int n_halfspaces,
                   const double *halfspaces,
                   const double *ray,
                   double *value)
{
    int rval = 0;

    *value = -INFINITY;

    for(int i = 0; i < n_halfspaces; i++)
    {
        const double *h = &halfspaces[2 * i];
        double v = ray[0] * h[0] + ray[1] * h[1];
        *value = DOUBLE_max(v, *value);
    }

CLEANUP:
    return rval;
}

int LIFTING_2D_optimize_continuous(int n_halfspaces,
                                   const double *halfspaces,
                                   double alpha2,
                                   double *alpha1,
                                   double *value)
{
    int rval = 0;

    *value = -INFINITY;

    for(int r = 0; r < n_halfspaces; r++)
    {
        for(int s = r+1; s < n_halfspaces; s++)
        {
            const double *ar = &halfspaces[r * 2];
            const double *as = &halfspaces[s * 2];

            double xr = ar[0];
            double xs = as[0];
            if(DOUBLE_eq(xr, xs)) continue;

            double br = - xs  / (xr - xs);
            double bs = xr / (xr - xs);
            if(!DOUBLE_geq(br, 0)) continue;
            if(!DOUBLE_geq(bs, 0)) continue;

            double yr = ar[1] * alpha2;
            double ys = as[1] * alpha2;

            double obj = yr * br + ys * bs;

            if(DOUBLE_geq(obj, *value))
            {
                *value = obj;

                if(DOUBLE_neq(xr, 0))
                    *alpha1 = - (yr - obj) / xr;
                else
                    *alpha1 = - (ys - obj) / xs;
            }
        }
    }

CLEANUP:
    return rval;
}

int LIFTING_2D_lift_fixed(int n_halfspaces,
                          const double *halfspaces,
                          const double *ray,
                          double k1,
                          double *opt)
{
    int rval = 0;

    double k0;
    double value;

    rval = LIFTING_2D_optimize_continuous(n_halfspaces, halfspaces,
            ray[1] + k1, &k0, &value);
    abort_if(rval, "LIFTING_2D_optimize_continuous failed");

    double delta = k0 - ray[0];

    double r_ceil[2] = { ray[0] + ceil(delta), ray[1] + k1 };
    double r_floor[2] = { ray[0] + floor(delta), ray[1] + k1 };

    log_debug("    delta=%.6lf value=%.6lf\n", delta, value);
    log_debug("    r_ceil=%12.6lf %12.6lf\n", r_ceil[0], r_ceil[1]);
    log_debug("    r_floor=%12.6lf %12.6lf\n", r_floor[0], r_floor[1]);

    double value_ceil, value_floor;

    rval = LIFTING_2D_psi(n_halfspaces, halfspaces, r_ceil, &value_ceil);
    abort_if(rval, "LIFTING_2D_psi failed");

    rval = LIFTING_2D_psi(n_halfspaces, halfspaces, r_floor, &value_floor);
    abort_if(rval, "LIFTING_2D_psi failed");

    *opt = min(value_ceil, value_floor);

CLEANUP:
    return rval;
}

int LIFTING_2D_naive(int n_halfspaces,
                     const double *halfspaces,
                     const double *ray,
                     const int *lb,
                     const int *ub,
                     double *value)
{
    int rval = 0;

    *value = INFINITY;
    int best_k0 = 0;
    int best_k1 = 0;

    int margin = 10;

    for(int k0 = lb[0] - margin; k0 <= ub[0] + margin; k0++)
    {
        for(int k1 = lb[1] - margin; k1 <= ub[1] + margin; k1++)
        {
            double q[2] = { ray[0] + k0, ray[1] + k1 };
            double value_q;

            rval = LIFTING_2D_psi(n_halfspaces, halfspaces, q, &value_q);
            abort_if(rval, "LIFTING_2D_ps failed");

            if(value_q < *value)
            {
                best_k0 = k0;
                best_k1 = k1;
                *value = value_q;
                log_debug("  k=%6d %6d value=%12.6lf\n", k0, k1, *value);
            }
        }
    }

//    log_debug("      best_k=(%d %d) value=%.6lf\n", best_k0, best_k1, *value);

CLEANUP:
    return rval;
}

int LIFTING_2D_bound(int n_halfspaces,
                     const double *halfspaces,
                     const double *ray,
                     double *value)
{
    int rval = 0;

    double eta_star, eta_plus, eta_minus;
    double xi_plus, xi_minus;
    double m_plus, m_minus;
    double ignored;
    int k1 = 1;

    rval = LIFTING_2D_lift_fixed(n_halfspaces, halfspaces, ray, 0, &eta_star);
    abort_if(rval, "LIFTING_2D_lift_fixed failed");

    rval = LIFTING_2D_optimize_continuous(n_halfspaces, halfspaces, 1, &ignored,
            &xi_plus);
    abort_if(rval, "LIFTING_2D_optimize_continuous failed");

    rval = LIFTING_2D_optimize_continuous(n_halfspaces, halfspaces, -1,
            &ignored, &xi_minus);
    abort_if(rval, "LIFTING_2D_optimize_continuous failed");

    m_plus = m_minus = INFINITY;

    log_debug("Level 0:\n");
    log_debug("  eta star  = %.6lf\n", eta_star);
    log_debug("  next slack plus  = %.6lf\n", k1 + ray[1] - m_plus);
    log_debug("  next slack minus = %.6lf\n", k1 - ray[1] - m_minus);
    log_debug("  xi plus = %.6lf\n", xi_plus);
    log_debug("  xi minus = %.6lf\n", xi_minus);

    while((k1 <= fabs(ray[1])) || (k1 + ray[1] <= m_plus) || (k1 - ray[1] <= m_minus))
    {
        log_debug("Level %d:\n", k1);

        double h_plus, h_minus;

        rval = LIFTING_2D_optimize_continuous(n_halfspaces, halfspaces, ray[1]
                + k1, &ignored, &h_plus);
        abort_if(rval, "LIFTING_2D_optimize_continuous failed");

        rval = LIFTING_2D_optimize_continuous(n_halfspaces, halfspaces, ray[1]
                - k1, &ignored, &h_minus);
        abort_if(rval, "LIFTING_2D_optimize_continuous failed");

        rval = LIFTING_2D_lift_fixed(n_halfspaces, halfspaces, ray, k1,
                &eta_plus);
        abort_if(rval, "LIFTING_2D_lift_fixed failed");

        rval = LIFTING_2D_lift_fixed(n_halfspaces, halfspaces, ray, -k1,
                &eta_minus);
        abort_if(rval, "LIFTING_2D_lift_fixed failed");

        eta_star = fmin(eta_star, fmin(eta_plus, eta_minus));

        m_plus = ceil(eta_star / xi_plus);
        m_minus = ceil(eta_star / xi_minus);

        log_debug("  eta plus  = %.6lf\n", eta_plus);
        log_debug("  eta minus = %.6lf\n", eta_minus);
        log_debug("  eta star  = %.6lf\n", eta_star);
        log_debug("    h minus = %.6lf\n", h_minus);
        log_debug("     h plus = %.6lf\n", h_plus);

        k1++;

        log_debug("  next slack plus  = %.6lf\n", k1 + ray[1] - m_plus);
        log_debug("  next slack minus = %.6lf\n", k1 - ray[1] - m_minus);
    }

    log_debug("Done\n");

    *value = eta_star;

CLEANUP:
    return rval;
}
