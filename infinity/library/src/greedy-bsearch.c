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
#include <stdlib.h>

#include <multirow/cg.h>
#include <multirow/double.h>
#include <multirow/geometry.h>
#include <multirow/lp.h>
#include <multirow/util.h>

#include <infinity/greedy-bsearch.h>

int create_sfree_mip(int nrows,
                     int nrays,
                     const double *f,
                     const double *rays,
                     const double *bounds,
                     double e,
                     struct LP *lp)
{
    int rval = 0;

    double rhs;
    char sense;

    int rmatbeg = 0;
    int* rmatind = 0;
    double *rmatval = 0;

    rmatind = (int *) malloc((nrows + nrays) * sizeof(int));
    rmatval = (double *) malloc((nrows + nrays) * sizeof(double));
    abort_if(!rmatind, "could not allocate rmatind");
    abort_if(!rmatval, "could not allocate rmatval");

    rval = LP_create(lp, "greedy");
    abort_if(rval, "LP_create failed");

    // create x (basic) variables
    for (int i = 0; i < nrows; i++)
    {
        rval = LP_new_col(lp, 0, -MILP_INFINITY, MILP_INFINITY, 'I');
        abort_if(rval, "LP_new_col failed");
    }

    // create s (non-basic) variables
    for (int i = 0; i < nrays; i++)
    {
        rval = LP_new_col(lp, 1.0, 0, MILP_INFINITY, 'C');
        abort_if(rval, "LP_new_col failed");
    }

    // add constraint \sum_{i=1}^m s_i \leq 1
    sense = 'L';
    rhs = 1.0;

    for (int i = 0; i < nrays; i++)
    {
        rmatind[i] = i + nrows;
        rmatval[i] = 1.0;
    }

    rval = LP_add_rows(lp, 1, nrays, &rhs, &sense, &rmatbeg, rmatind, rmatval);
    abort_if(rval, "LP_add_rows failed");

    // add constraints x_i - \sum_{j=1}^m min{e,e_j} s_j R_ji = f_i
    for (int i = 0; i < nrows; i++)
    {
        int k = 0;
        sense = 'E';
        rhs = f[i];

        rmatind[k] = i;
        rmatval[k] = 1.0;
        k++;

        for (int j = 0; j < nrays; j++)
        {
            rmatind[k] = j + nrows;
            rmatval[k] = -rays[nrows * j + i] * fmin(e, bounds[j]);
            k++;
        }

        rval = LP_add_rows(lp, 1, nrays + 1, &rhs, &sense, &rmatbeg, rmatind,
                           rmatval);
        abort_if(rval, "LP_add_rows failed");
    }

CLEANUP:
    if (rmatind) free(rmatind);
    if (rmatval) free(rmatval);
    return rval;
}

int GREEDY_BSEARCH_compute_bounds(int nrows,
                                  int nrays,
                                  const double *f,
                                  const double *rays,
                                  double *bounds)
{
    int rval = 0;

    struct LP lp;
    double e_upper = 2 * GREEDY_BIG_E;
    double e_lower = 0.0;

    int cplex_count = 0;
    double cplex_time = 0;

    int iteration_count = 0;

    double *x = 0;

    x = (double *) malloc((nrays + nrows) * sizeof(double));
    abort_if(!x, "could not allocate x");

    for (int i = 0; i < nrays; i++)
        bounds[i] = GREEDY_BIG_E;

    for (int it = 0;; it++)
    {
        abort_if(it > 2*nrays, "stuck in an infinite loop");

        log_verbose("Starting iteration %d...\n", it);

        iteration_count++;


        int solution_found = 0;
        int inner_count = 0;

        while (fabs(e_upper - e_lower) > GREEDY_MAX_GAP)
        {
            inner_count++;

            double e = (e_upper + e_lower) / 2;
            log_verbose("    e=%.12lf\n", e);

            rval = LP_open(&lp);
            abort_if(rval, "LP_open failed");

            rval = create_sfree_mip(nrows, nrays, f, rays, bounds, e, &lp);
            abort_if(rval, "create_sfree_mip failed");

            if_verbose_level
            {
                rval = LP_write(&lp, "greedy.lp");
                abort_if(rval, "LP_write failed");
            }

            int infeasible;
            cplex_count++;

            double initial_time = get_user_time();

            log_verbose("    Optimizing...\n");
            rval = LP_optimize(&lp, &infeasible);
            if (rval)
            {
                // Workaround for CPLEX bug. If CPLEX tell us that this problem
                // is unbounded, we disable presolve and try again.
                LP_free(&lp);
                LP_open(&lp);

                rval = create_sfree_mip(nrows, nrays, f, rays, bounds, e, &lp);
                abort_if(rval, "create_sfree_mip failed");

                LP_disable_presolve(&lp);

                rval = LP_optimize(&lp, &infeasible);
                abort_if(rval, "LP_optimize failed");
            }

            cplex_time += get_user_time() - initial_time;

            if (infeasible)
            {
                e_lower = e;
                log_verbose("    infeasible\n");
                if (e > GREEDY_BIG_E-1)
                {
                    LP_free(&lp);
                    goto OUT;
                }
            }
            else
            {
                log_verbose("    feasible\n");
                e_upper = e;
                solution_found = 1;

                rval = LP_get_x(&lp, x);
                abort_if(rval, "LP_get_x failed");
            }

            LP_free(&lp);
        }

        if (solution_found)
        {
            for (int j = 0; j < nrays; j++)
            {
                if (!DOUBLE_geq(x[nrows + j], 0.001)) continue;
                bounds[j] = fmin(bounds[j] * 0.99, e_lower * 0.99);
            }
        }

        log_verbose("    %d iterations  %12.8lf gap\n", inner_count, e_upper -
                    e_lower);

        e_lower = e_upper;
        e_upper = 2 * GREEDY_BIG_E;
    }

OUT:
    log_debug("    %6d IPs (%.2lfms per call, %.2lfs total)\n", cplex_count,
                cplex_time * 1000.0 / cplex_count, cplex_time);

    for(int i = 0; i < nrays; i++)
        abort_if(DOUBLE_iszero(bounds[i]), "bounds should be positive");

    if_verbose_level
    {
        time_printf("Bounds:\n");
        for (int k = 0; k < nrays; k++)
            time_printf("    %12.8lf  %12.8lf\n", k, bounds[k], 1 / bounds[k]);
    }

CLEANUP:
    if (x) free(x);
    return rval;
}
