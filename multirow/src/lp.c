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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ilcplex/cplex.h>

#include <multirow/lp.h>
#include <multirow/util.h>
#include <multirow/double.h>

int LP_open(struct LP *lp)
{
    int rval = 0;

    lp->cplex_lp = 0;
    lp->cplex_env = CPXopenCPLEX(&rval);
    abort_if(rval, "CPXopenCPLEX failed");

    CPXsetintparam(lp->cplex_env, CPX_PARAM_DATACHECK, CPX_ON);
    CPXsetintparam(lp->cplex_env, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
    CPXsetlogfile(lp->cplex_env, 0);

CLEANUP:
    return rval;
}

void LP_disable_presolve(struct LP *lp)
{
    CPXsetintparam(lp->cplex_env, CPX_PARAM_PREIND, CPX_OFF);
    CPXsetintparam(lp->cplex_env, CPX_PARAM_REPEATPRESOLVE, CPX_OFF);
}

int LP_create(struct LP *lp,
              const char *name)
{
    int rval = 0;
    char nambuf[MAX_NAME_LENGTH];

    abort_if(!lp->cplex_env, "cplex_env is null");

    strncpy(nambuf, name, MAX_NAME_LENGTH);
    nambuf[MAX_NAME_LENGTH - 1] = '\0';

    lp->cplex_lp = CPXcreateprob(lp->cplex_env, &rval, nambuf);
    abort_if(rval, "CPXcreateprob failed");

CLEANUP:
    return rval;
}

void LP_free(struct LP *lp)
{
    if (!lp) return;
    if (!lp->cplex_env) return;

    if (lp->cplex_lp)
        CPXfreeprob(lp->cplex_env, &(lp->cplex_lp));

    CPXcloseCPLEX(&lp->cplex_env);
    lp->cplex_env = 0;
}

void LP_free_row(struct Row *row)
{
    if (!row) return;
    if (row->pi) free(row->pi);
    if (row->indices) free(row->indices);
}

int LP_add_rows(struct LP *lp,
                int newrows,
                int newnz,
                double *rhs,
                char *sense,
                int *rmatbeg,
                int *rmatind,
                double *rmatval)
{
    int rval = 0;

    rval = CPXaddrows(lp->cplex_env, lp->cplex_lp, 0, newrows, newnz, rhs,
                      sense, rmatbeg, rmatind, rmatval, 0, 0);
    abort_if(rval, "CPXaddrows failed");

CLEANUP:
    return rval;
}

int LP_add_row(struct LP *lp,
               struct Row *row)
{
    int rval = 0;

    char sense = 'L';
    int rmatbeg = 0;

    rval = CPXaddrows(lp->cplex_env, lp->cplex_lp, 0, 1, row->nz, &row->pi_zero,
                      &sense, &rmatbeg, row->indices, row->pi, 0, 0);
    abort_if(rval, "CPXaddrows failed");

CLEANUP:
    return rval;
}

int LP_add_cols(struct LP *lp,
                int newcols,
                int newnz,
                double *obj,
                int *cmatbeg,
                int *cmatind,
                double *cmatval,
                double *lb,
                double *ub)
{
    int rval = 0;

    rval = CPXaddcols(lp->cplex_env, lp->cplex_lp, newcols, newnz, obj, cmatbeg,
                      cmatind, cmatval, lb, ub, 0);
    abort_if(rval, "CPXaddcols failed");

CLEANUP:
    return rval;
}

int LP_change_bound(struct LP *lp,
                    int col,
                    char lower_or_upper,
                    double bnd)
{
    int rval = 0;

    rval = CPXchgbds(lp->cplex_env, lp->cplex_lp, 1, &col, &lower_or_upper,
                     &bnd);
    abort_if(rval, "CPXchgbds failed");

CLEANUP:
    return rval;
}

int LP_get_obj_val(struct LP *lp,
                   double *obj)
{
    int rval = 0;

    rval = CPXgetobjval(lp->cplex_env, lp->cplex_lp, obj);
    abort_if(rval, "CPXgetobjval failed");

CLEANUP:
    return rval;
}

int LP_get_x(struct LP *lp,
             double *x)
{
    int rval = 0;

    int ncols = CPXgetnumcols(lp->cplex_env, lp->cplex_lp);
    abort_if(!ncols, "No columns in LP");

    rval = CPXgetx(lp->cplex_env, lp->cplex_lp, x, 0, ncols - 1);
    abort_if(rval, "CPXgetx failed");

CLEANUP:
    return rval;
}

int LP_get_y(struct LP *lp,
             double *y)
{
    int rval = 0;

    int nrows = CPXgetnumrows(lp->cplex_env, lp->cplex_lp);
    abort_if(!nrows, "No rows in LP");

    rval = CPXgetpi(lp->cplex_env, lp->cplex_lp, y, 0, nrows - 1);
    abort_iff(rval, "CPXgetpi failed (errno = %d)", rval);

CLEANUP:
    return rval;
}

int LP_get_num_cols(const struct LP *lp)
{
    return CPXgetnumcols(lp->cplex_env, lp->cplex_lp);
}

int LP_get_num_rows(const struct LP *lp)
{
    return CPXgetnumrows(lp->cplex_env, lp->cplex_lp);
}

int LP_new_row(struct LP *lp,
               char sense,
               double rhs)
{
    int rval = 0;

    rval = CPXnewrows(lp->cplex_env, lp->cplex_lp, 1, &rhs, &sense, 0, 0);
    abort_if(rval, "CPXnewrows failed");

CLEANUP:
    return rval;
}

int LP_new_col(struct LP *lp,
               double obj,
               double lb,
               double ub,
               char type)
{
    int rval = 0;

    rval = CPXnewcols(lp->cplex_env, lp->cplex_lp, 1, &obj, &lb, &ub, &type, 0);
    abort_if(rval, "CPXnewcols failed");

CLEANUP:
    return rval;
}

int LP_optimize(struct LP *lp,
                int *infeasible)
{
    int rval = 0, solstat;

    *infeasible = 0;

//    LP_SOLVE_COUNT++;

    if_verbose_level
    {
        int numrows = CPXgetnumrows(lp->cplex_env, lp->cplex_lp);
        int numcols = CPXgetnumcols(lp->cplex_env, lp->cplex_lp);

        time_printf("Optimizing MILP (%d rows %d cols)...\n", numrows, numcols);
    }

    double initial_time = get_user_time();
    int problem_type = CPXgetprobtype(lp->cplex_env, lp->cplex_lp);

    switch (problem_type)
    {
    case CPXPROB_LP:
        rval = CPXdualopt(lp->cplex_env, lp->cplex_lp);
        log_verbose("  dual opt\n");
        abort_if(rval, "CPXdualopt failed");
        break;

    case CPXPROB_MILP:
    case CPXPROB_FIXEDMILP:
        rval = CPXmipopt(lp->cplex_env, lp->cplex_lp);
        log_verbose("  mip opt\n");
        abort_if(rval, "CPXmipopt failed");
        break;

    default:
        abort_iff(1, "Invalid problem type: %d", problem_type);
    }

    solstat = CPXgetstat(lp->cplex_env, lp->cplex_lp);
    switch (solstat)
    {
    case CPX_STAT_INFEASIBLE:
    case CPXMIP_INFEASIBLE:
        log_verbose("    infeasible\n");
        *infeasible = 1;
        goto CLEANUP;

    case CPXMIP_OPTIMAL:
    case CPXMIP_OPTIMAL_INFEAS:
    case CPXMIP_OPTIMAL_TOL:
    case CPX_STAT_OPTIMAL:
    case CPX_STAT_OPTIMAL_INFEAS:
    case CPX_STAT_UNBOUNDED:
        break;

    default:
        abort_iff(1, "Invalid solution status: %d", solstat);
    }

    double objval;
    rval = LP_get_obj_val(lp, &objval);
    abort_if(rval, "LP_get_obj_val failed");

    log_verbose("    obj val = %.4lf\n", objval);
    log_verbose("    time = %.4lf\n", get_user_time() - initial_time);

CLEANUP:
    return rval;
}

int LP_write(struct LP *lp,
             const char *fname)
{
    int rval = 0;
    char nambuf[MAX_NAME_LENGTH];

    FILE *f = fopen(fname, "w");
    abort_iff(!f, "could not open file %s", fname);
    fclose(f);

    strncpy(nambuf, fname, MAX_NAME_LENGTH);
    nambuf[MAX_NAME_LENGTH - 1] = '\0';

    log_info("Writing LP to file %s...\n", fname);
    rval = CPXwriteprob(lp->cplex_env, lp->cplex_lp, nambuf, "RLP");
    abort_if(rval, "CPXwriteprob failed");

CLEANUP:
    return rval;
}

int LP_read_problem(struct LP *lp,
                    const char *filename)
{
    int rval = 0;

    log_info("Reading problem %s...\n", filename);

    rval = CPXreadcopyprob(lp->cplex_env, lp->cplex_lp, filename, 0);
    abort_if(rval, "CPXreadcopyprob failed");

CLEANUP:
    return rval;
}

int LP_relax(struct LP *lp)
{
    int rval = 0;

    rval = CPXchgprobtype(lp->cplex_env, lp->cplex_lp, CPXPROB_LP);
    abort_if(rval, "CPXchgprobtype failed");

CLEANUP:
    return rval;
}

int LP_get_column_types(struct LP *lp,
                        char *column_types)
{
    int rval = 0;
    char *cplex_ctype = 0;

    int ncols = LP_get_num_cols(lp);

    cplex_ctype = (char *) malloc(ncols * sizeof(char));
    abort_if(!cplex_ctype, "could not allocate cplex_ctype");

    rval = CPXgetctype(lp->cplex_env, lp->cplex_lp, cplex_ctype, 0, ncols - 1);
    abort_if(rval, "CPXgetctype failed");

    for (int i = 0; i < ncols; i++)
    {
        switch (cplex_ctype[i])
        {
        case CPX_BINARY:
        case CPX_INTEGER:
            column_types[i] = MILP_INTEGER;
            break;

        case CPX_CONTINUOUS:
            column_types[i] = MILP_CONTINUOUS;
            break;

        default:
            abort_iff(1, "Invalid column type: %d", cplex_ctype[i]);
        }
    }

CLEANUP:
    if (cplex_ctype) free(cplex_ctype);
    return rval;
}

int LP_get_tableau(struct LP *lp,
                   struct Row **rows,
                   int *cstat,
                   int *rstat,
                   double *ub,
                   double *lb)
{
    int rval = 0;

    int nrows = LP_get_num_rows(lp);
    int ncols = LP_get_num_cols(lp);

    int *head = 0;
    double *rhs = 0;
    double *pi = 0;
    double *slacks = 0;

    rhs = (double *) malloc(nrows * sizeof(double));
    head = (int *) malloc(nrows * sizeof(int));
    slacks = (double*) malloc(nrows * sizeof(double));

    abort_if(!head, "could not allocate head");
    abort_if(!rhs, "could not allocate rhs");
    abort_if(!slacks, "could not allocate slacks");

    rval = CPXgetbhead(lp->cplex_env, lp->cplex_lp, head, rhs);
    abort_if(rval, "CPXgetbhead failed");

    rval = CPXgetslack(lp->cplex_env, lp->cplex_lp, slacks, 0, nrows-1);
    abort_if(rval, "CPXgetslack failed");

    rval = CPXgetub(lp->cplex_env, lp->cplex_lp, ub, 0, ncols - 1);
    abort_if(rval, "CPXgetub failed");

    rval = CPXgetlb(lp->cplex_env, lp->cplex_lp, lb, 0, ncols - 1);
    abort_if(rval, "CPXgetlb failed");

    rval = CPXgetbase(lp->cplex_env, lp->cplex_lp, cstat, rstat);
    abort_if(rval, "CPXgetbase failed");

    for (int i = 0; i < nrows; i++)
        abort_if(!DOUBLE_iszero(slacks[i]), "model contains slack variables");

    pi = (double *) malloc(ncols * sizeof(double));
    abort_if(!pi, "could not allocate pi");

    for (int i = 0; i < nrows; i++)
    {
        int nz = 0;
        rows[i] = (struct Row *) malloc(sizeof(struct Row));
        abort_if(!rows[i], "could not allocate rows[i]");

        struct Row *row = rows[i];

        rval = CPXbinvarow(lp->cplex_env, lp->cplex_lp, i, pi);
        abort_if(rval, "CPXbinvarow failed");

        for (int j = 0; j < ncols; j++)
        {
            if (fabs(pi[j]) < EPSILON)
                continue;

            if (cstat[j] == CPX_AT_LOWER)
                rhs[i] += lb[j] * pi[j];

            if (cstat[j] == CPX_AT_UPPER)
                rhs[i] += ub[j] * pi[j];

            nz++;
        }

        row->nz = nz;
        row->pi_zero = rhs[i];
        row->head = head[i];

        row->pi = (double *) malloc(nz * sizeof(double));
        row->indices = (int *) malloc(nz * sizeof(int));

        abort_if(!row->pi, "could not allocate row->pi");
        abort_if(!row->indices, "could not allocate row->indices");

        if (fabs(row->pi_zero) < EPSILON)
            row->pi_zero = 0;

        int k = 0;
        for (int j = 0; j < ncols; j++)
        {
            if (fabs(pi[j]) < EPSILON)
                continue;

            row->pi[k] = pi[j];
            row->indices[k++] = j;
        }

        rval = LP_flip_row_coefficients(cstat, lb, ub, row);
        abort_if(rval, "LP_flip_row_coefficients failed");
    }

CLEANUP:
    if(slacks) free(slacks);
    if (head) free(head);
    if (rhs) free(rhs);
    if (pi) free(pi);
    return rval;
}

int LP_flip_row_coefficients(int *cstat,
                             double *lb,
                             double *ub,
                             struct Row *cut)
{
    for (int j = 0; j < cut->nz; j++)
    {
        int idx = cut->indices[j];
        double pij = cut->pi[j];

        if (cstat[idx] == CPX_AT_LOWER)
            cut->pi_zero -= lb[idx] * pij;

        if (cstat[idx] == CPX_AT_UPPER)
        {
            cut->pi_zero -= ub[idx] * pij;
            pij = -pij;
        }

        cut->indices[j] = idx;
        cut->pi[j] = pij;
    }

    return 0;
}

int LP_unflip_row_coefficients(int *cstat,
                               double *lb,
                               double *ub,
                               struct Row *cut)
{
    for (int j = 0; j < cut->nz; j++)
    {
        int idx = cut->indices[j];
        double pij = cut->pi[j];

        if (cstat[idx] == CPX_AT_LOWER)
            cut->pi_zero += lb[idx] * pij;

        if (cstat[idx] == CPX_AT_UPPER)
        {
            pij = -pij;
            cut->pi_zero += ub[idx] * pij;
        }

        cut->indices[j] = idx;
        cut->pi[j] = pij;
    }

    return 0;
}

void LP_print_row(struct Row *row)
{
    for (int i = 0; i < row->nz; i++)
        printf("%.2lfx%d ", row->pi[i], row->indices[i]);

    printf(" <= %.2lf\n", row->pi_zero);
}

int LP_read_solution(struct LP *lp,
                     char *filename,
                     double *x)
{
    int rval = 0;
    int ncols = LP_get_num_cols(lp);

    log_info("Reading solution from %s\n", filename);

    FILE *fsol = fopen(filename, "rb");
    abort_iff(!fsol, "Could not open file %s", filename)

    for (int i = 0; i < ncols; i++)
    {
        int count = fscanf(fsol, "%le", &x[i]);
        abort_iff(count != 1, "Unexpected EOF when reading %s", filename);
    }

    fclose(fsol);

CLEANUP:
    return rval;
}

int LP_write_solution(struct LP *lp,
                      char *filename)
{
    int rval = 0;
    long ncols = LP_get_num_cols(lp);

    double *x = 0;

    x = (double *) malloc(ncols * sizeof(double));
    abort_if(!x, "could not allocate x");

    rval = LP_get_x(lp, x);
    abort_if(rval, "LP_get_x failed");

    log_info("Writing solution to file %s\n", filename);

    FILE *fsol = fopen(filename, "w");
    abort_iff(!fsol, "Could not open file %s", filename)

    for (int i = 0; i < ncols; i++)
        fprintf(fsol, "%.12e\n", x[i]);

    fclose(fsol);

CLEANUP:
    if (x) free(x);
    return rval;
}

int LP_write_basis(struct LP *lp,
                   char *filename)
{
    int rval = 0;

    log_info("Writing basis to file %s...\n", filename);
    rval = CPXmbasewrite(lp->cplex_env, lp->cplex_lp, filename);
    abort_if(rval, "CPXmbasewrite failed");

CLEANUP:
    return rval;
}

int LP_read_basis(struct LP *lp,
                  char *filename)
{
    int rval = 0;

    log_info("Reading basis from file %s...\n", filename);
    rval = CPXreadcopybase(lp->cplex_env, lp->cplex_lp, filename);
    abort_if(rval, "CPXreadcopybase failed");

CLEANUP:
    return rval;
}


int LP_write_sage_file(struct LP *lp,
                       double *x,
                       char *filename)
{
    int rval = 0;

    FILE *fsage = 0;


    int ncols = LP_get_num_cols(lp);
    int nrows = LP_get_num_rows(lp);

    int nz;
    int dummy;

    int *indices = 0;
    double *rhs = 0;
    double *values = 0;
    int *head = 0;

    fsage = fopen(filename, "w");
    abort_if(!fsage, "fopen failed");

    indices = (int*) malloc(ncols * sizeof(int));
    values = (double*) malloc(ncols * sizeof(double));
    rhs = (double*) malloc(nrows * sizeof(double));
    head = (int*) malloc(nrows * sizeof(int));

    abort_if(!indices, "could not allocate indices");
    abort_if(!values, "could not allocate values");
    abort_if(!rhs, "could not allocate rhs");
    abort_if(!head, "could not allocate head");

    rval = CPXgetrhs(lp->cplex_env, lp->cplex_lp, rhs, 0, nrows-1);
    abort_if(rval, "CPXgetrhs failed");

    fprintf(fsage, "b = vector([");
    for (int i = 0; i < nrows; i++)
        fprintf(fsage, "  %.8lf,\n", rhs[i]);
    fprintf(fsage, "])\n");

    fprintf(fsage, "A = matrix(QQ, %d, %d, sparse=True)\n", nrows, ncols);

    for (int i = 0; i < nrows; i++)
    {
        rval = CPXgetrows(lp->cplex_env, lp->cplex_lp, &nz, &dummy, indices,
                          values, ncols, &dummy, i, i);
        abort_if(rval, "CPXgetrows failed");

        for (int k = 0; k < nz; k++)
            fprintf(fsage, "A[%5d,%5d] = %20.8lf\n", i, indices[k], values[k]);
    }

    rval = CPXgetbhead(lp->cplex_env, lp->cplex_lp, head, rhs);
    abort_if(rval, "CPXgetbhead failed");

    fprintf(fsage, "B=[");
    for (int i = 0; i < nrows; i++)
        fprintf(fsage, " %d,\n", head[i]);
    fprintf(fsage, "]\n");

    fprintf(fsage, "x=vector([");
    for (int i = 0; i < ncols; i++)
        fprintf(fsage, "  %.8lf,\n", x[i]);
    fprintf(fsage, "])\n");

    fclose(fsage);

    CLEANUP:
    if(indices) free(indices);
    if(values) free(values);
    if(rhs) free(rhs);
    if(head) free(head);
    return rval;
}

int LP_change_rhs(struct LP *lp, int index, double value)
{
    int rval = 0;

    rval = CPXchgrhs(lp->cplex_env, lp->cplex_lp, 1, &index, &value);
    abort_if(rval, "CPXchgrhs failed");

CLEANUP:
    return rval;
}

int LP_init_row(struct Row *row, int nz_capacity)
{
    int rval = 0;

    row->nz = 0;
    row->head = 0;
    row->pi_zero = 0;

    row->pi = (double *) malloc(nz_capacity * sizeof(double));
    row->indices = (int *) malloc(nz_capacity * sizeof(int));
    abort_if(!row->pi, "could not allocate row->pi");
    abort_if(!row->indices, "could not allocate row->indices");

    CLEANUP:
    return rval;
}