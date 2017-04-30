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

#ifndef _PROJECT_LP_H_
#define _PROJECT_LP_H_

#include <ilcplex/cplex.h>
#include "params.h"

#define MILP_INTEGER 0
#define MILP_CONTINUOUS 1
#define MILP_INFINITY CPX_INFBOUND

struct LP
{
    CPXENVptr cplex_env;
    CPXLPptr cplex_lp;
};

struct Row
{
    int nz;
    int head;
    double pi_zero;
    double *pi;
    int *indices;
};

static const int MAX_NAME_LENGTH = 100;

int LP_open(struct LP *lp);

int LP_create(struct LP *lp,
              const char *name);

int LP_write(struct LP *lp,
             const char *fname);

void LP_free(struct LP *lp);

void LP_free_row(struct Row *row);

int LP_new_row(struct LP *lp,
               char sense,
               double rhs);

int LP_new_col(struct LP *lp,
               double obj,
               double lb,
               double ub,
               char type);

int LP_add_row(struct LP *lp,
               struct Row *row);

int LP_add_rows(struct LP *lp,
                int newrows,
                int newnz,
                double *rhs,
                char *sense,
                int *rmatbeg,
                int *rmatind,
                double *rmatval);

int LP_add_cols(struct LP *lp,
                int newcols,
                int newnz,
                double *obj,
                int *cmatbeg,
                int *cmatind,
                double *cmatval,
                double *lb,
                double *ub);

int LP_change_bound(struct LP *lp,
                    int col,
                    char lower_or_upper,
                    double bnd);

int LP_optimize(struct LP *lp,
                int *infeasible);

int LP_get_obj_val(struct LP *lp,
                   double *obj);

int LP_get_x(struct LP *lp,
             double *x);

int LP_get_y(struct LP *lp,
             double *y);

int LP_get_num_cols(const struct LP *lp);

int LP_get_num_rows(const struct LP *lp);

int LP_read_problem(struct LP *lp,
                    const char *filename);

int LP_relax(struct LP *lp);

int LP_get_column_types(struct LP *lp,
                        char *column_types);

int LP_get_tableau(struct LP *lp,
                   struct Row **rows,
                   int *cstat,
                   int *rstat,
                   double *ub,
                   double *lb);

int LP_unflip_row_coefficients(int *cstat,
                               double *lb,
                               double *ub,
                               struct Row *cut);

void LP_print_row(struct Row *row);

int LP_flip_row_coefficients(int *cstat,
                             double *lb,
                             double *ub,
                             struct Row *cut);

int LP_read_solution(struct LP *lp,
                     char *filename,
                     double *x);

void LP_disable_presolve(struct LP *lp);

int LP_write_solution(struct LP *lp,
                      char *filename);

int LP_write_basis(struct LP *lp,
                   char *filename);

int LP_read_basis(struct LP *lp,
                  char *filename);

int LP_write_sage_file(struct LP *lp,
                       double *x,
                       char *filename);

int LP_change_rhs(struct LP *lp, int index, double value);

int LP_init_row(struct Row *row, int nz_capacity);

#endif
