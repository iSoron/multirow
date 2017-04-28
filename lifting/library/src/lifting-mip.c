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

#include <multirow/util.h>
#include <multirow/lp.h>
#include <lifting/lifting-mip.h>

double MIP_TIME_OPTIMIZE = 0;
double MIP_TIME_CREATE = 0;

/**
 * Returns non-zero if the cplex_env part of struct LP is initialized.
 *
 * @param lp  pointer to a struct LP
 */
int LP_is_open(struct LP *lp)
{
    return (lp->cplex_env != 0);
}

/**
 * Frees the cplex_lp part of struct LP, leaving the environment open.
 *
 * @param lp  pointer to a struct LP
 */
void LP_destroy(struct LP *lp)
{
    if(!lp) return;
    if(!lp->cplex_env) return;

    if(lp->cplex_lp)
        CPXfreeprob(lp->cplex_env, &(lp->cplex_lp));

    lp->cplex_lp = 0;
}

/**
 * Static data for LIFTING_2D_mip().
 *
 * Warning: Because of this, LIFTING_2D_mip() is not thread safe.
 */
static struct LP LIFTING_2D_mip_lp = {0, 0};

/**
 * Opens the cplex environment of the struct LP used by LIFTING_2D_mip().
 *
 * This is called automatically by LIFTING_2D_mip().
 */
int LIFTING_2D_mip_init()
{
    struct LP *lp = &LIFTING_2D_mip_lp;

    if(!LP_is_open(lp))
    {
        if(LP_open(lp)) return (-1);
    }

    return (0);
}

/**
 * Closes the cplex environment of the struct LP used by LIFTING_2D_mip().
 *
 * This function should be called by the user of LIFTING_2D_mip(), if
 * cleaning up the environment is desired.
 */
void LIFTING_2D_mip_cleanup()
{
    struct LP *lp = &LIFTING_2D_mip_lp;

    LP_free(lp);
}

/*
 * Computes the lifting coefficient of a ray by formulating the lifting
 * problem as a MIP.
 *
 * @param[in]  n_halfspaces  number of facets of the lattice-free set used
 * @param[in]  halfspaces    description of the lattice-free set used
 * @param[in]  ray           ray to lift
 * @param[out] value         lifting coefficient
 */
int LIFTING_2D_mip(int n_halfspaces,
                   const double *halfspaces,
                   const double *ray,
                   double *value)
{
    struct LP *lp = &LIFTING_2D_mip_lp;
    int rval = 0;

    double initial_time = get_user_time();

    rval = LIFTING_2D_mip_init();
    abort_if(rval, "LIFTING_2D_mip_init failed");

    rval = LP_create(lp, "lifting2d");
    abort_if(rval, "LP_create failed");

    rval = LP_new_col(lp, 1.0, 0.0, MILP_INFINITY, 'C');
    rval |= LP_new_col(lp, 0.0, -MILP_INFINITY, MILP_INFINITY, 'I');
    rval |= LP_new_col(lp, 0.0, -MILP_INFINITY, MILP_INFINITY, 'I');
    abort_if(rval, "LP_new_col failed");

    double rhs;
    char sense = 'G';
    int beg = 0;
    int ind[3] = {0, 1, 2};
    double val[3];

    val[0] = 1.0;

    for(int i = 0; i < n_halfspaces; i++)
    {
        double a0 = halfspaces[i * 2 + 0];
        double a1 = halfspaces[i * 2 + 1];

        rhs = a0 * ray[0] + a1 * ray[1];

        val[1] = -a0;
        val[2] = -a1;

        rval |= LP_add_rows(lp, 1, 3, &rhs, &sense, &beg, ind, val);
    }

    MIP_TIME_CREATE += get_user_time() - initial_time;

    abort_if(rval, "LP_add_rows failed");

    int infeasible;

    initial_time = get_user_time();

    rval = LP_optimize(lp, &infeasible);
    abort_if(rval, "LP_optimize failed");
    abort_if(infeasible, "LIFTING_2D_mip infeasible");

    MIP_TIME_OPTIMIZE += get_user_time() - initial_time;

    double obj;

    rval = LP_get_obj_val(lp, &obj);
    abort_if(rval, "LP_get_obj_val failed");

    *value = obj;

CLEANUP:
    LP_destroy(lp);
    return (rval);
}
