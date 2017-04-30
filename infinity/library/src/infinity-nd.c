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
#include <multirow/util.h>

#include <infinity/infinity-nd.h>

static long lp_count = 0;
static double lp_time = 0;

static long epsilon_lp_count = 0;
static double epsilon_lp_time = 0;

static long tight_lp_count = 0;
static double tight_lp_time = 0;

static long violated_lp_count = 0;
static double violated_lp_time = 0;

static long sfree_mip_count = 0;
static double sfree_mip_time = 0;

static long scale_ahull_lp_count = 0;
static double scale_ahull_lp_time = 0;

static int create_sfree_mip(int nrows,
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

static int create_find_epsilon_lp(int nrows,
                                  int nrays,
                                  const double *f,
                                  const double *rays,
                                  const int *t,
                                  const int *rx,
                                  const double *x,
                                  const double *beta,
                                  struct LP *lp)
{
    int rval = 0;

    double rhs;
    char sense;

    int nz = 0;
    int *map = 0;
    int rmatbeg = 0;
    int *rmatind = 0;
    double *rmatval = 0;
    int rx_count = 0;

    map = (int *) malloc(nrays * sizeof(int));
    abort_if(!map, "could not allocate map");

    rmatind = (int *) malloc((nrays + 1 + nrows) * sizeof(int));
    rmatval = (double *) malloc((nrays + 1 + nrows) * sizeof(double));
    abort_if(!rmatind, "could not allocate rmatind");
    abort_if(!rmatval, "could not allocate rmatval");

    rval = LP_create(lp, "find_epsilon");
    abort_if(rval, "LP_create failed");

    // create lambda variables
    for(int i = 0; i < nrays + 1; i++)
    {
        if(i < nrays && !rx[i]) continue;

        double pi = 0.0;
        double lb = -MILP_INFINITY;

        if(i < nrays && !t[i])
        {
            pi = 1.0;
            lb = 0.0;
        }

        rval = LP_new_col(lp, pi, lb, MILP_INFINITY, 'C');
        abort_if(rval, "LP_new_col failed");

        if(i < nrays)
            map[i] = rx_count++;
    }

    log_verbose("rx_count=%d\n", rx_count);

    // create y variables
    for(int i = 0; i < nrows; i++)
    {
        rval = LP_new_col(lp, 0.0, -MILP_INFINITY, MILP_INFINITY, 'C');
        abort_if(rval, "LP_new_col failed");
    }

    // create constraint y = \lambda_x x + \sum_{t \in T} \lambda_r (f + \beta_r r)
    for(int j = 0; j < nrows; j++)
    {
        sense = 'E';
        rhs = 0.0;
        nz = 0;

        for(int i = 0; i < nrays; i++)
        {
            if(!t[i]) continue;
            const double *ri = &rays[i * nrows];

            rmatind[nz] = map[i];
            rmatval[nz] = f[j] + beta[i] * ri[j];
            nz++;
        }

        rmatind[nz] = rx_count;
        rmatval[nz] = x[j];
        nz++;

        rmatind[nz] = rx_count + j + 1;
        rmatval[nz] = -1.0;
        nz++;

        rval = LP_add_rows(lp, 1, nz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
        abort_if(rval, "LP_add_rows failed");
    }

    // create constraint y = f + \sum_{r \in Rx \setminus T) \lambda_r r
    for(int j = 0; j < nrows; j++)
    {
        sense = 'E';
        rhs = f[j];
        nz = 0;

        for(int i = 0; i < nrays; i++)
        {
            if(!rx[i] || t[i]) continue;
            const double *ri = &rays[i * nrows];

            rmatind[nz] = map[i];
            rmatval[nz] = -ri[j];
            nz++;
        }

        rmatind[nz] = rx_count + j + 1;
        rmatval[nz] = 1.0;
        nz++;

        rval = LP_add_rows(lp, 1, nz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
        abort_if(rval, "LP_add_rows failed");
    }

    // create constraint \sum_{r \in T} \lambda_r + \lambda_x = 1
    sense = 'E';
    rhs = 1.0;
    nz = 0;

    for(int i = 0; i < nrays; i++)
    {
        if(!t[i]) continue;
        rmatind[nz] = map[i];
        rmatval[nz] = 1.0;
        nz++;
    }

    rmatind[nz] = rx_count;
    rmatval[nz] = 1.0;
    nz++;

    rval = LP_add_rows(lp, 1, nz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
    abort_if(rval, "LP_add_rows failed");

    rval = LP_relax(lp);
    abort_if(rval, "LP_relax failed");

    //rval = LP_write(lp, "find-epsilon.lp");
    //abort_if(rval, "LP_write failed");


CLEANUP:
    if(map) free(map);
    if(rmatind) free(rmatind);
    if(rmatval) free(rmatval);
    return rval;
}

static int create_tight_rays_lp(int nrows,
                                int nrays,
                                const double *f,
                                const double *rays,
                                const double *x,
                                const double *beta,
                                double epsilon,
                                double delta,
                                struct LP *lp)
{
    int rval = 0;

    double rhs;
    char sense;

    int rmatbeg = 0;
    int *rmatind = 0;
    double *rmatval = 0;

    rmatind = (int *) malloc(nrays * sizeof(int));
    rmatval = (double *) malloc(nrays * sizeof(double));
    abort_if(!rmatind, "could not allocate rmatind");
    abort_if(!rmatval, "could not allocate rmatval");

    rval = LP_create(lp, "tight_rays");
    abort_if(rval, "LP_create failed");

    // create lambda variables
    for(int i = 0; i < nrays; i++)
    {
        rval = LP_new_col(lp, 1.0, 0.0, MILP_INFINITY, 'C');
        abort_if(rval, "LP_new_col failed");
    }

    // create s variables
    for(int i = 0; i < nrays; i++)
    {
        rval = LP_new_col(lp, 0.0, 0.0, MILP_INFINITY, 'C');
        abort_if(rval, "LP_new_col failed");
    }

    // create constraint x = f + \sum_{r \in R} min{e, beta[r]} * r * s_r
    for(int j = 0; j < nrows; j++)
    {
        sense = 'E';
        rhs = x[j] - f[j];

        for(int i = 0; i < nrays; i++)
        {
            const double *ri = &rays[i * nrows];
            rmatind[i] = nrays + i;
            rmatval[i] = min(epsilon, beta[i]) * ri[j];
            if(DOUBLE_iszero(rmatval[i])) rmatval[i] = 0.0;
        }

        rval = LP_add_rows(lp, 1, nrays, &rhs, &sense, &rmatbeg, rmatind,
                rmatval);
        abort_if(rval, "LP_add_rows failed");
    }

    // create constraint \sum_{r \in R} s_r = 1
    sense = 'E';
    rhs = 1.0;

    for(int i = 0; i < nrays; i++)
    {
        rmatind[i] = nrays + i;
        rmatval[i] = 1.0;
    }

    rval = LP_add_rows(lp, 1, nrays, &rhs, &sense, &rmatbeg, rmatind, rmatval);
    abort_if(rval, "LP_add_rows failed");

    // create constraints \lambda_r + s_r \geq \delta
    for(int i = 0; i < nrays; i++)
    {
        sense = 'G';
        rhs = delta;

        rmatind[0] = i;
        rmatval[0] = 1.0;

        rmatind[1] = nrays + i;
        rmatval[1] = 1.0;

        rval = LP_add_rows(lp, 1, 2, &rhs, &sense, &rmatbeg, rmatind, rmatval);
        abort_if(rval, "LP_add_rows failed");
    }

    rval = LP_relax(lp);
    abort_if(rval, "LP_relax failed");

    //rval = LP_write(lp, "tight-rays.lp");
    //abort_if(rval, "LP_write failed");


CLEANUP:
    if(rmatind) free(rmatind);
    if(rmatval) free(rmatval);
    return rval;
}

static int create_violated_cone_lp(int nrows,
                                   int nrays,
                                   const double *f,
                                   const double *rays,
                                   const double *x,
                                   const double *beta,
                                   double epsilon,
                                   struct LP *lp)
{
    int rval = 0;

    double rhs;
    char sense;

    int rmatbeg = 0;
    int *rmatind = 0;
    double *rmatval = 0;

    rmatind = (int *) malloc(nrays * sizeof(int));
    rmatval = (double *) malloc(nrays * sizeof(double));
    abort_if(!rmatind, "could not allocate rmatind");
    abort_if(!rmatval, "could not allocate rmatval");

    rval = LP_create(lp, "violated_cone");
    abort_if(rval, "LP_create failed");

    // create s variables
    for(int i = 0; i < nrays; i++)
    {
        rval = LP_new_col(lp, 1.0, 0.0, MILP_INFINITY, 'C');
        abort_if(rval, "LP_new_col failed");
    }

    // create constraint x = f + \sum(min{e, beta[r]} * r * s_r)
    for(int j = 0; j < nrows; j++)
    {
        sense = 'E';
        rhs = x[j] - f[j];

        for(int i = 0; i < nrays; i++)
        {
            const double *ri = &rays[i * nrows];
            rmatind[i] = i;
            rmatval[i] = min(epsilon, beta[i]) * ri[j];
            if(DOUBLE_iszero(rmatval[i])) rmatval[i] = 0.0;
        }

        rval = LP_add_rows(lp, 1, nrays, &rhs, &sense, &rmatbeg, rmatind,
                rmatval);
        abort_if(rval, "LP_add_rows failed");
    }

    rval = LP_relax(lp);
    abort_if(rval, "LP_relax failed");

    //rval = LP_write(lp, "violated-cone.lp");
    //abort_if(rval, "LP_write failed");
    //UTIL_pause();

CLEANUP:
    if(rmatind) free(rmatind);
    if(rmatval) free(rmatval);
    return rval;
}

static int create_scale_to_ahull_lp(int nrows,
                                    int nrays,
                                    const double *rays,
                                    const int *rx,
                                    const double *beta,
                                    double epsilon,
                                    const double *d,
                                    struct LP *lp)
{
    int rval = 0;

    double rhs;
    char sense;
    int nz;

    int rmatbeg = 0;
    int *rmatind = 0;
    double *rmatval = 0;

    rmatind = (int *) malloc((nrays + 1) * sizeof(int));
    rmatval = (double *) malloc((nrays + 1) * sizeof(double));
    abort_if(!rmatind, "could not allocate rmatind");
    abort_if(!rmatval, "could not allocate rmatval");

    rval = LP_create(lp, "scale_to_ahull");
    abort_if(rval, "LP_create failed");

    // create alpha variable
    rval = LP_new_col(lp, 1.0, 0.0, MILP_INFINITY, 'C');
    abort_if(rval, "LP_new_col failed");


    // create lambda variables
    for(int i = 0; i < nrays; i++)
    {
        rval = LP_new_col(lp, 0.0, -MILP_INFINITY, MILP_INFINITY, 'C');
        abort_if(rval, "LP_new_col failed");
    }

    // create constraint \sum_{r \in R_x} min(e, beta[r]) * r * \lambda_r = \alpha * d
    for(int j = 0; j < nrows; j++)
    {
        sense = 'E';
        rhs = 0.0;
        nz = 0;

        rmatind[nz] = 0;
        rmatval[nz] = d[j];
        nz++;

        for(int i = 0; i < nrays; i++)
        {
            if(!rx[i]) continue;

            const double *ri = &rays[i * nrows];
            rmatind[nz] = 1 + i;
            rmatval[nz] = -min(epsilon, beta[i]) * ri[j];
            if(DOUBLE_iszero(rmatval[nz])) continue;

            nz++;
        }

        rval = LP_add_rows(lp, 1, nz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
        abort_if(rval, "LP_add_rows failed");
    }

    // create constraint \sum_{r \in R_x} \lambda_r = 1
    sense = 'E';
    rhs = 1.0;
    nz = 0;

    for(int i = 0; i < nrays; i++)
    {
        if(!rx[i]) continue;

        rmatind[nz] = 1 + i;
        rmatval[nz] = 1.0;
        nz++;
    }

    rval = LP_add_rows(lp, 1, nz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
    abort_if(rval, "LP_add_rows failed");

    rval = LP_relax(lp);
    abort_if(rval, "LP_relax failed");

    //rval = LP_write(lp, "scale-to-ahull.lp");
    //abort_if(rval, "LP_write failed");
    //UTIL_pause();

CLEANUP:
    if(rmatind) free(rmatind);
    if(rmatval) free(rmatval);
    return rval;
}

static int find_interior_point_enum(const int nrows,
                                    const int nrays,
                                    const double *f,
                                    const double *rays,
                                    const double *beta,
                                    const double epsilon,
                                    double *x,
                                    int *found)
{
    int rval = 0;

    int M = 1;
    struct LP lp;
    double best_value = INFINITY;
    double *beta2 = 0;

    abort_if(nrows != 3, "not implemented");

    beta2 = (double *) malloc(nrays * sizeof(double));
    abort_if(!beta2, "could not allocate beta2");

    for(int i = 0; i < nrays; i++)
    {
        beta2[i] = fmin(epsilon, beta[i]);
    }

    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    struct ConvLFreeSet lfree;
    lfree.f = (double*) f;
    lfree.nrows = nrows;
    lfree.rays.dim = nrows;
    lfree.rays.nrays = nrays;
    lfree.rays.values = (double*) rays;
    lfree.beta = beta2;

    rval = INFINITY_create_psi_lp(&lfree, &lp);
    abort_if(rval, "INFINITY_create_psi_lp failed");

    *found = 0;

    for(int x1 = -M; x1 <= M; x1++)
        for(int x2 = -M; x2 <= M; x2++)
            for(int x3 = -M; x3 <= M; x3++)
            {
                double value;
                double q[3] = {x1 - f[0],
                               x2 - f[1],
                               x3 - f[2]};

                rval = INFINITY_psi(nrows, q, 1, &lp, &value);
                abort_if(rval, "INFINITY_psi failed");

                if(value < best_value)
                {
                    best_value = value;
                    x[0] = x1;
                    x[1] = x2;
                    x[2] = x3;
                }
            }

    if(best_value < 0.999) *found = 1;

CLEANUP:
    if(beta2) free(beta2);
    LP_free(&lp);
    return rval;
}

static int find_interior_point_cplex(const int nrows,
                                     const int nrays,
                                     const double *f,
                                     const double *rays,
                                     const double *beta,
                                     const double epsilon,
                                     double *x,
                                     int *found)
{
    int rval = 0;
    struct LP lp;
    double initial_time;
    int infeasible;
    double objval;

    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    lp_count++;
    sfree_mip_count++;
    initial_time = get_user_time();

    rval = create_sfree_mip(nrows, nrays, f, rays, beta, epsilon, &lp);
    abort_if(rval, "greate_sfree_mip failed");

    log_debug("  solving sfree mip...\n");
    rval = LP_optimize(&lp, &infeasible);
    if(rval == ERR_MIP_TIMEOUT) goto CLEANUP;
    abort_if(rval, "LP_optimize failed");

    if(infeasible)
    {
        log_debug("  mip is infeasible. Stopping.\n");
        *found = 0;
        goto CLEANUP;
    }

    rval = LP_get_x(&lp, x);
    abort_if(rval, "LP_get_x failed");

    if_verbose_level
    {
        for(int i = 0; i < nrows; i++)
                log_verbose("    x%d = %.8lf\n", i, x[i]);

        for(int i = 0; i < nrays; i++)
            if(x[i + nrows] > 0.00001)
                    log_verbose("    t%d = %.8lf\n", i, x[i + nrows]);
    }

    rval = LP_get_obj_val(&lp, &objval);
    abort_if(rval, "LP_get_obj_val failed");

    log_debug("    obj = %.8lf\n", objval);

    if(objval >= 0.999)
    {
        log_debug("  set is lattice-free\n");
        *found = 0;
        goto CLEANUP;
    }

    *found = 1;

    sfree_mip_time += get_user_time() - initial_time;
    lp_time += get_user_time() - initial_time;

CLEANUP:
    LP_free(&lp);
    return rval;
}

static int cone_bound(int nrows,
                      int nrays,
                      const double *f,
                      const double *rays,
                      const int *rx,
                      const double *x,
                      const double *beta,
                      double *epsilon)
{
    log_verbose("      find_epsilon\n");
    int rval = 0;

    int *t = 0;
    long it = 0;

    t = (int *) malloc(nrays * sizeof(int));
    abort_if(!t, "could not allocate t");

    for(int i = 0; i < nrays; i++)
        t[i] = 0;

    while(1)
    {
        it++;
        log_verbose("Starting iteration %d...\n", it);

        struct LP lp;

        double initial_time = get_user_time();
        rval = LP_open(&lp);
        abort_if(rval, "LP_open failed");

        rval = create_find_epsilon_lp(nrows, nrays, f, rays, t, rx, x, beta,
                &lp);
        abort_if(rval, "create_find_epsilon_lp failed");

        int infeasible;
        rval = LP_optimize(&lp, &infeasible);
        abort_if(rval, "LP_optimize failed");

        lp_count++;
        lp_time += get_user_time() - initial_time;

        epsilon_lp_count++;
        epsilon_lp_time += get_user_time() - initial_time;

        if(infeasible)
        {
            *epsilon = INFINITY;
            log_verbose("  infeasible\n");
            LP_free(&lp);
            goto CLEANUP;
        }

        double obj;
        rval = LP_get_obj_val(&lp, &obj);
        abort_if(rval, "LP_get_obj_val failed");

        log_verbose("    obj=%.6lf\n", obj);

        for(int i = 0; i < nrays; i++)
        {
            if(!rx[i]) continue;
            log_verbose("  beta[%d]=%.6lf\n", i, beta[i]);
        }

        LP_free(&lp);

        double e_min = INFINITY;
        double e_max = -INFINITY;

        for(int i = 0; i < nrays; i++)
        {
            if(!rx[i] || t[i]) continue;
            e_min = min(e_min, beta[i]);
        }

        for(int i = 0; i < nrays; i++)
        {
            if(!rx[i]) continue;
            e_max = fmax(e_max, beta[i]);
        }

        log_verbose("  e_max=%.6lf\n", e_max);
        log_verbose("  e_min=%.6lf\n", e_min);

        if(DOUBLE_leq(obj, e_min))
        {
            if(DOUBLE_geq(obj, e_max))
                *epsilon = INFINITY;
            else
                *epsilon = obj;

            goto CLEANUP;
        }
        else
        {
            for(int i = 0; i < nrays; i++)
                if(rx[i] && DOUBLE_eq(beta[i], e_min))
                    t[i] = 1;
        }
    }

CLEANUP:
    log_verbose("  e=%.6lf\n", *epsilon);
    if(t) free(t);
    return rval;
}

static int find_tight_rays(int nrows,
                           int nrays,
                           const double *f,
                           const double *rays,
                           const double *x,
                           const double *beta,
                           double epsilon,
                           int *tx)
{
    log_verbose("      find_tight_rays\n");
    int rval = 0;
    int infeasible = 0;
    double initial_time = 0;
    const double delta = 0.001;

    struct LP lp;
    double *sbar = 0;

    sbar = (double *) malloc(2 * nrays * sizeof(double));
    abort_if(!sbar, "could not allocate sbar");

    initial_time = get_user_time();
    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    rval = create_tight_rays_lp(nrows, nrays, f, rays, x, beta, epsilon, delta,
            &lp);
    abort_if(rval, "create_tight_rays_lp failed");

    rval = LP_optimize(&lp, &infeasible);
    abort_if(rval, "LP_optimize failed");

    lp_count++;
    lp_time += get_user_time() - initial_time;

    tight_lp_count++;
    tight_lp_time += get_user_time() - initial_time;

    abort_if(infeasible, "tight_rays_lp is infeasible");

    rval = LP_get_x(&lp, sbar);
    abort_if(rval, "LP_get_x failed");

    for(int i = 0; i < nrays; i++)
        tx[i] = DOUBLE_iszero(sbar[i]);

    for(int i = 0; i < nrays; i++)
            log_verbose("  tx[%d]=%d\n", i, tx[i]);

CLEANUP:
    if(sbar) free(sbar);
    LP_free(&lp);
    return rval;
}

static int find_violated_cone(int nrows,
                              int nrays,
                              const double *f,
                              const double *rays,
                              const double *x,
                              const double *beta,
                              double epsilon,
                              int *rx,
                              double *sbar,
                              int *violated_found)
{
    log_verbose("      find_violated_cone\n");
    int rval = 0;

    struct LP lp;

    double initial_time = get_user_time();

    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    rval = create_violated_cone_lp(nrows, nrays, f, rays, x, beta, epsilon,
            &lp);
    abort_if(rval, "create_violated_cone_lp failed");

    int infeasible;
    rval = LP_optimize(&lp, &infeasible);
    abort_if(rval, "LP_optimize failed");

    lp_count++;
    lp_time += get_user_time() - initial_time;

    violated_lp_count++;
    violated_lp_time += get_user_time() - initial_time;

    rval = LP_get_x(&lp, sbar);
    abort_if(rval, "LP_get_x failed");

    for(int i = 0; i < nrays; i++)
        rx[i] = 0;

    if(infeasible)
        goto CLEANUP;

    double obj;
    rval = LP_get_obj_val(&lp, &obj);
    abort_if(rval, "LP_get_obj_val failed");

    log_verbose("  o=%.8lf\n", obj);

    if(DOUBLE_geq(obj, 0.999))
    {
        *violated_found = 0;
    }
    else
    {
        *violated_found = 1;

        log_verbose("Violated cone found\n");
        log_verbose("  f=%.8lf %.8lf\n", f[0], f[1]);
        log_verbose("  x=%.8lf %.8lf\n", x[0], x[1]);

        for(int i = 0; i < nrays; i++)
        {
            rx[i] = (sbar[i] > 1e-9);

            if(rx[i]) if_verbose_level
            {
                double m = min(epsilon, beta[i]);
                const double *r = &rays[nrows * i];
                time_printf("  r[%d]=%.8lf %.8lf\n", i, r[0], r[1]);
                time_printf("  r[%d]=%.8lf %.8lf\n", i, m * r[0], m * r[1]);
            }
        }
    }

CLEANUP:
    LP_free(&lp);
    return rval;
}

static int bound(int nrows,
                 int nrays,
                 const double *f,
                 const double *rays,
                 const double *x,
                 const double *beta,
                 double *epsilon,
                 int *tx)
{
    int rval = 0;

    int found;
    int *rx = 0;
    double *fbar = 0;
    double *sbar = 0;

    double prev_epsilon;
    int count = 0;
    *epsilon = GREEDY_BIG_E;

    rx = (int *) malloc(nrays * sizeof(int));
    fbar = (double *) malloc(nrows * sizeof(double));
    sbar = (double *) malloc(nrays * sizeof(double));

    abort_if(!rx, "could not allocate rx");
    abort_if(!fbar, "could not allocate fbar");
    abort_if(!sbar, "could not allocate sbar");

    while(1)
    {
        count++;
        abort_if(count > 100, "infinite loop");

        rval = find_violated_cone(nrows, nrays, f, rays, x, beta, *epsilon, rx,
                sbar, &found);
        abort_if(rval, "find_violated_cone failed");

        if(!found) break;

        for(int i = 0; i < nrows; i++)
            fbar[i] = x[i];

        for(int j = 0; j < nrays; j++)
        {
            if(!rx[j]) continue;
            const double *r = &rays[nrows * j];

            for(int i = 0; i < nrows; i++)
                fbar[i] -= min(*epsilon, beta[j]) * r[i] * sbar[j];
        }

        log_verbose("%.12lf %.12lf\n", f[0], f[1]);
        log_verbose("%.12lf %.12lf\n", fbar[0], fbar[1]);

        prev_epsilon = *epsilon;

        rval = cone_bound(nrows, nrays, fbar, rays, rx, x, beta, epsilon);
        abort_if(rval, "cone_bound failed");

        log_verbose("        e=%.12lf\n", *epsilon);
        abort_if(prev_epsilon < *epsilon, "epsilon should never increase");
    }

    for(int i = 0; i < nrays; i++)
        tx[i] = 0;

    if(DOUBLE_geq(*epsilon, GREEDY_BIG_E))
    {
        *epsilon = INFINITY;
        goto CLEANUP;
    }
    else
    {
        rval = find_tight_rays(nrows, nrays, fbar, rays, x, beta, *epsilon, tx);
        abort_if(rval, "find_tight_rays failed");
    }

CLEANUP:
    if(sbar) free(sbar);
    if(fbar) free(fbar);
    if(rx) free(rx);
    return rval;
}

static int scale_to_ahull(int nrows,
                          int nrays,
                          const double *rays,
                          const int *rx,
                          const double *beta,
                          double epsilon,
                          const double *d,
                          double *alpha)
{
    log_verbose("      scale_to_ahull\n");
    int rval = 0;
    *alpha = INFINITY;

    struct LP lp;
    double *x = 0;
    double initial_time;

    x = (double *) malloc((nrays + 1) * sizeof(double));
    abort_if(!x, "could not allocate x");

    initial_time = get_user_time();
    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    rval = create_scale_to_ahull_lp(nrows, nrays, rays, rx, beta, epsilon, d,
            &lp);
    abort_if(rval, "create_scale_to_ahull_lp failed");

    int infeasible;
    rval = LP_optimize(&lp, &infeasible);
    abort_if(rval, "LP_optimize failed");

    lp_count++;
    lp_time += get_user_time() - initial_time;

    scale_ahull_lp_count++;
    scale_ahull_lp_time += get_user_time() - initial_time;

    if(infeasible)
        goto CLEANUP;

    rval = LP_get_x(&lp, x);
    abort_if(rval, "LP_get_x failed");

    *alpha = x[0];

CLEANUP:
    if(x) free(x);
    LP_free(&lp);
    return rval;
}

#ifndef TEST_SOURCE

int INFINITY_create_psi_lp(const struct ConvLFreeSet *lfree, struct LP *lp)
{
    int rval = 0;

    int rmatbeg = 0;
    int *rmatind = 0;
    double *rmatval = 0;

    int nrows = lfree->nrows;
    int nrays = lfree->rays.nrays;
    double *rays = lfree->rays.values;
    double *beta = lfree->beta;

    rmatind = (int *) malloc(nrays * sizeof(int));
    rmatval = (double *) malloc(nrays * sizeof(double));
    abort_if(!rmatind, "could not allocate rmatind");
    abort_if(!rmatval, "could not allocate rmatval");

    rval = LP_create(lp, "psi");
    abort_if(rval, "LP_create failed");

    // create lambda variables
    for(int i = 0; i < nrays; i++)
    {
        rval = LP_new_col(lp, 1.0, 0, MILP_INFINITY, 'C');
        abort_if(rval, "LP_new_col failed");
    }

    // create constraint 0 = \sum_{i=1}^m \lambda_i r^i_j beta_i
    for(int j = 0; j < nrows; j++)
    {
        char sense = 'E';
        double rhs = 0;
        int nz = 0;

        for(int i = 0; i < nrays; i++)
        {
            const double *ri = &rays[i * nrows];
            rmatind[nz] = i;
            rmatval[nz] = ri[j] * beta[i];
            nz++;
        }

        rval = LP_add_rows(lp, 1, nz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
        abort_if(rval, "LP_add_rows failed");
    }

    rval = LP_relax(lp);
    abort_if(rval, "LP_relax failed");

    if_verbose_level
    {
        rval = LP_write(lp, "psi.lp");
        abort_if(rval, "LP_write failed");
    }

CLEANUP:
    if(rmatind) free(rmatind);
    if(rmatval) free(rmatval);
    return rval;
}

int INFINITY_pi(const int nrows,
                const double *q,
                const double q_scale,
                struct LP *lp,
                double *value)
{
    int rval = 0;

    abort_if(nrows > 3, "not implemented");

    double best_value = 1e100;

    if(nrows == 2)
    {
        int M = 2;
        for(int k0 = -M; k0 <= M; k0++)
            for(int k1 = -M; k1 <= M; k1++)
            {
                double v;
                double qk[] = {frac(q[0] * q_scale) + k0,
                               frac(q[1] * q_scale) + k1};

                rval = INFINITY_psi(nrows, qk, 1, lp, &v);
                abort_if(rval, "INFINITY_psi failed");

                best_value = min(best_value, v);
            }
    }

    if(nrows == 3)
    {
        int M = 2;
        for(int k0 = -M; k0 <= M; k0++)
            for(int k1 = -M; k1 <= M; k1++)
                for(int k2 = -M; k2 <= M; k2++)
                {
                    double v;
                    double qk[] = {frac(q[0] * q_scale) + k0,
                                   frac(q[1] * q_scale) + k1,
                                   frac(q[2] * q_scale) + k2};

                    rval = INFINITY_psi(nrows, qk, 1, lp, &v);
                    abort_if(rval, "INFINITY_psi failed");

                    best_value = min(best_value, v);
                }
    }

    *value = best_value;

CLEANUP:
    return rval;
}

/*
 * Given a point f, a list of rays r1,...,rm, some real numbers b1,...,bm, a
 * point q, and a real number q_scale, this function evaluates psi_B(q *
 * q_scale), where B is the convex hull of {f + ri * bi}_i=1^m.
 */
int INFINITY_psi(const int nrows,
                 const double *q,
                 const double q_scale,
                 struct LP *lp,
                 double *value)
{
    int rval = 0;

    for(int j = 0; j < nrows; j++)
    {
        rval = LP_change_rhs(lp, j, q[j] * q_scale);
        abort_if(rval, "LP_change_rhs failed");
    }

    int infeasible;
    rval = LP_optimize(lp, &infeasible);
    abort_if(rval, "LP_optimize failed");

    if(infeasible)
    {
        *value = INFINITY;
    }
    else
    {
        rval = LP_get_obj_val(lp, value);
        abort_if(rval, "LP_get_obj_val failed");
    }

CLEANUP:
    return rval;
}

int INFINITY_ND_generate_lfree(const struct MultiRowModel *model,
                               struct ConvLFreeSet *lfree)
{
    int rval = 0;
    int nrows = model->nrows;
    int nrays = model->rays.nrays;

    double *f = lfree->f;
    double *rays = lfree->rays.values;
    double *beta = lfree->beta;

    lfree->nrows = model->nrows;
    lfree->rays.nrays = nrays;
    memcpy(rays, model->rays.values, nrays * nrows * sizeof(double));
    memcpy(f, model->f, nrows * sizeof(double));

    double *x = 0;

    int *t = 0;
    int *tx = 0;

    t = (int *) malloc(nrays * sizeof(int));
    tx = (int *) malloc(nrays * sizeof(int));
    abort_if(!t, "could not allocate t");
    abort_if(!tx, "could not allocate tx");

    x = (double *) malloc((nrows + nrays) * sizeof(double));
    abort_if(!x, "could not allocate x");

    for(int i = 0; i < nrays; i++)
        beta[i] = GREEDY_BIG_E;

    int it = 0;

    //lp_time = 0;
    //lp_count = 0;

    //epsilon_lp_count = 0;
    //epsilon_lp_time = 0;
    //
    //sfree_mip_count = 0;
    //sfree_mip_time = 0;

    //tight_lp_count = 0;
    //tight_lp_time = 0;

    //violated_lp_count = 0;
    //violated_lp_time = 0;

    //scale_ahull_lp_time = 0;
    //scale_ahull_lp_count = 0;

    long x_count = 0;
    double epsilon;

    for(int i = 0; i < nrows; i++)
            log_verbose("  f[%d] = %.12lf\n", i, f[i]);

    while(1)
    {
        it++;
        abort_if(it > 10 * nrays, "infinite loop");

        log_debug("Starting iteration %d...\n", it);
        epsilon = INFINITY;

        for(int i = 0; i < nrays; i++)
            t[i] = 0;

        for(int i = 0; i < nrows; i++)
            x[i] = 0;

        while(1)
        {
            log_debug("  epsilon = %.8lf\n", epsilon);

            int found = 0;

            if(nrows == 3)
            {
                rval = find_interior_point_enum(nrows, nrays, f,
                        model->rays.values, beta, epsilon, x, &found);
                abort_if(rval, "find_interior_point_enum failed");
            }

            if(!found)
            {
                rval = find_interior_point_cplex(nrows, nrays, f,
                        model->rays.values, beta, epsilon, x, &found);
                if(rval == ERR_MIP_TIMEOUT) goto CLEANUP;
                abort_if(rval, "find_interior_point_cplex failed");
                if(!found) break;
            }

            log_debug("    found interior point:\n");
            for(int i = 0; i < nrows; i++) log_debug("        %.2lf\n", x[i]);

            x_count++;
            abort_if(x_count > 1000, "infinite loop");

            double epsilon_x;
            rval = bound(nrows, nrays, f, model->rays.values, x, beta,
                    &epsilon_x, tx);
            abort_if(rval, "bound failed");
//            epsilon_x *= 0.999;

            if(isinf(epsilon_x)) break;

            log_debug("    epsilon_x = %.8lf\n", epsilon_x);

            if(DOUBLE_eq(epsilon_x, epsilon))
            {
                for(int i = 0; i < nrays; i++)
                    if(tx[i]) t[i] = 1;
            }
            else if(epsilon_x < epsilon)
            {
                epsilon = epsilon_x;
                for(int i = 0; i < nrays; i++)
                    t[i] = tx[i];
            }
        }

        if(isinf(epsilon))
            break;

        int skip_ahull = 0;

        for(int i = 0; i < nrays; i++)
        {
            if(t[i])
            {
                beta[i] = min(beta[i], epsilon);
//                beta[i] *= 0.999;
            }
            else if(!skip_ahull)
            {
                double alpha;
                const double *d = LFREE_get_ray(&model->rays, i);

                rval = scale_to_ahull(nrows, nrays, model->rays.values, t, beta,
                        epsilon, d, &alpha);
                abort_if(rval, "scale_to_ahull failed");

                if(DOUBLE_iszero(alpha))
                {
                    skip_ahull = 1;
                    continue;
                }

                beta[i] = min(beta[i], alpha);
//                beta[i] *= 0.999;
            }

            log_debug("  beta[%2d] = %.4lf\n", i, beta[i]);
        }

        log_debug("epsilon = %.6lf\n", epsilon);
    }

    log_debug("    %6ld lattice points, %ld iterations\n", x_count, it);

    log_debug("    %6ld MIPs (%.2lf ms per call, %.0lf ms total)\n", lp_count,
            lp_time * 1000 / lp_count, lp_time * 1000);

    log_debug(
            "          %6ld S-free MIPs (%.2lf ms per call, %.0lf ms total)\n",
            sfree_mip_count, sfree_mip_time * 1000 / sfree_mip_count,
            sfree_mip_time * 1000);

    log_debug(
            "          %6ld epsilon LPs (%.2lf ms per call, %.0lf ms total)\n",
            epsilon_lp_count, epsilon_lp_time * 1000 / epsilon_lp_count,
            epsilon_lp_time * 1000);

    log_debug(
            "          %6ld tight-rays LPs (%.2lf ms per call, %.0lf ms total)\n",
            tight_lp_count, tight_lp_time * 1000 / tight_lp_count,
            tight_lp_time * 1000);

    log_debug(
            "          %6ld violated-cone LPs (%.2lf ms per call, %.0lf ms total)\n",
            violated_lp_count, violated_lp_time * 1000 / violated_lp_count,
            violated_lp_time * 1000);

    log_debug(
            "          %6ld scale-to-ahull LPs (%.2lf ms per call, %.0lf ms total)\n",
            scale_ahull_lp_count,
            scale_ahull_lp_time * 1000 / scale_ahull_lp_count,
            scale_ahull_lp_time * 1000);

CLEANUP:
    if(x) free(x);
    if(t) free(t);
    if(tx) free(tx);
    return rval;
}

#endif // TEST_SOURCE