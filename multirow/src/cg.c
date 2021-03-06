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

#include <stdlib.h>
#include <math.h>

#include <multirow/params.h>
#include <multirow/cg.h>
#include <multirow/double.h>
#include <multirow/stats.h>
#include <multirow/util.h>

static double cg_initial_time;

static int select_rows(struct CG *cg, int *row_selected)
{
    int rval = 0;

    long histogram[100] = {0};

    for (int i = 0; i < cg->nrows; i++)
    {
        row_selected[i] = 1;
        struct Row *r1 = cg->tableau_rows[i];

        if (r1->head < 0 || cg->column_types[r1->head] != MILP_INTEGER)
        {
            row_selected[i] = 0;
            continue;
        }

        double dist = fabs(0.5 - frac(r1->pi_zero));
        for (int k = 0; k < 100; k++)
            if (dist <= (k / 100.0))
                histogram[k]++;
    }

    double dist_cutoff = 0.5;
    for (int k = 0; k < 100; k++)
    {
        if (histogram[k] > MAX_SELECTED_ROWS)
        {
            dist_cutoff = k / 100.0;
            break;
        }
    }

    int selected_count = 0;
    for (int i = 0; i < cg->nrows; i++)
    {
        if (!row_selected[i]) continue;

        struct Row *r1 = cg->tableau_rows[i];
        double dist = fabs(0.5 - frac(r1->pi_zero));

        if (dist > dist_cutoff || selected_count >= MAX_SELECTED_ROWS)
        {
            row_selected[i] = 0;
            continue;
        }

        selected_count++;
    }

CLEANUP:
    return rval;
}

static int
next_combination(int n, int k, int *array, int inc_index, int *finished)
{
    int i;
    int rval = 0;

    array[inc_index]++;
    *finished = 0;

    for (i = inc_index; i < k; i++)
    {
        if (array[i] < n - i)
            break;

        if (i + 1 < k)
            array[i + 1]++;
        else
        {
            *finished = 1;
            return 0;
        }
    }

    for (int j = 1; j <= i; j++)
        array[i - j] = array[i] + j;

    return 0;
}

static int ray_norm(int dim, const double *ray, double *norm)
{
    *norm = 0;

    for (int i = 0; i < dim; i++)
        *norm += fabs(ray[i]);

    return 0;
}

/*
 * Checks whether two rays are parallel. Also returns the scaling factor,
 * in case they are.
 */
static int check_rays_parallel(int dim,
                               const double *r1,
                               const double *r2,
                               int *match,
                               double *scale)
{
    int rval = 0;

    *match = 1;

    double r1_norm, r2_norm;
    int iszero_r1, iszero_r2;

    rval = ray_norm(dim, r1, &r1_norm);
    abort_if(rval, "ray_norm failed");

    rval = ray_norm(dim, r2, &r2_norm);
    abort_if(rval, "ray_norm failed");

    iszero_r1 = DOUBLE_iszero(r1_norm);
    iszero_r2 = DOUBLE_iszero(r2_norm);

    if (iszero_r1 || iszero_r2)
    {
        *match = iszero_r1 && iszero_r2;
    }
    else
    {
        for (int i = 0; i < dim; i++)
        {
            if (DOUBLE_sgn(r1[i]) != DOUBLE_sgn(r2[i]) ||
                    DOUBLE_neq(r1[i] * r2_norm, r2[i] * r1_norm))
            {
                *match = 0;
                break;
            }

            log_verbose("        %.12lf equals %.12lf\n", r1[i] * r2_norm,
                    r2[i] * r1_norm);
        }
    }

    (*scale) = r1_norm / r2_norm;

CLEANUP:
    return rval;
}

static int
evaluate_row_pair(const struct Row *row1, const struct Row *row2, double *score)
{
    int rval = 0;

    int i1 = 0;
    int i2 = 0;
    int hit = 0;

    while (i1 < row1->nz || i2 < row2->nz)
    {
        int idx1 = INT_MAX;
        int idx2 = INT_MAX;

        if (i1 < row1->nz)
            idx1 = row1->indices[i1];

        if (i2 < row2->nz)
            idx2 = row2->indices[i2];

        if (idx1 == idx2) hit++;

        int idx_min = min(idx1, idx2);

        if (idx1 == idx_min)
            i1++;

        if (idx2 == idx_min)
            i2++;
    }

    *score = (hit * 2.0) / (row1->nz + row2->nz);

CLEANUP:
    return rval;
}

static int copy_solution(struct CG *cg, double *from, double **to)
{
    int rval = 0;
    int ncols = LP_get_num_cols(cg->lp);

    if (!*to)
    {
        *to = (double *) malloc(ncols * sizeof(double));
        abort_if(!*to, "could not allocate cg->integral_solution");
    }

    memcpy(*to, from, sizeof(double) * ncols);

CLEANUP:
    return rval;
}

static int check_cut(struct CG *cg, struct Row *cut)
{
    int rval = 0;

    if_verbose_level
    {
        time_printf("Checking cut:\n");
        for (int i = 0; i < cut->nz; i++)
        {
            if (DOUBLE_iszero(cg->integral_solution[cut->indices[i]])) continue;
            time_printf("    %12.8lf x%d", cut->pi[i], cut->indices[i]);
            if (cg->integral_solution)
                printf(" (=%12.8lf)", cg->integral_solution[cut->indices[i]]);
            printf("\n");
        }
        time_printf("     <= %.4lf\n", cut->pi_zero);
    }

    if (cg->basic_solution)
    {
        log_verbose("Basic solution check:\n");

        double lhs = CG_replace_x(cut, cg->basic_solution);

        log_verbose("    %.8lf > %.8lf\n", lhs, cut->pi_zero);
        abort_iff(!DOUBLE_geq(lhs, cut->pi_zero), "Cut fails to cut "
                "basic solution: %12.8lf < %12.8lf", lhs, cut->pi_zero);
    }

    if (cg->integral_solution)
    {
        log_verbose("Integral solution check:\n");

        double lhs = CG_replace_x(cut, cg->integral_solution);

        log_verbose("    %.8lf <= %.8lf\n", lhs, cut->pi_zero);
        abort_iff(!DOUBLE_leq(lhs, cut->pi_zero), "Cut cuts off known integral "
                "solution: %12.8lf >= %12.8lf", lhs, cut->pi_zero);
    }

CLEANUP:
    return rval;
}

static int add_cut(struct CG *cg, struct Row *cut, int *ignored)
{
    int rval;
    double *x = 0;
    double lhs;

    rval = LP_unflip_row_coefficients(cg->cstat, cg->lb, cg->ub, cut);
    abort_if(rval, "LP_unflip_row_coefficients failed");

    rval = check_cut(cg, cut);
    abort_if(rval, "check_cut failed");

    lhs = CG_replace_x(cut, cg->current_solution);

    *ignored = 0;
    STATS_increment_generated_cuts();

    if (DOUBLE_leq(lhs, cut->pi_zero))
    {
        log_verbose("Ignoring cut (%12.8lf <= %12.8lf)\n", lhs, cut->pi_zero);
        *ignored = 1;
    }
    else
    {
        rval = LP_add_row(cg->lp, cut);
        abort_if(rval, "LP_add_row failed");

        int infeasible;
        log_verbose("Reoptimizing...\n");
        rval = LP_optimize(cg->lp, &infeasible);
        abort_if(rval, "LP_optimize failed");

        double obj;
        rval = LP_get_obj_val(cg->lp, &obj);
        abort_if(rval, "LP_get_obj_val failed");

        if (DOUBLE_neq(obj, cg->last_obj_value))
            log_info("    opt = %lf\n", obj);

        cg->last_obj_value = obj;

        x = (double *) malloc(cg->ncols * sizeof(double));
        abort_if(!x, "could not allocate x");

        rval = LP_get_x(cg->lp, x);
        abort_if(rval, "LP_get_x failed");

        rval = copy_solution(cg, x, &cg->current_solution);
        abort_if(rval, "copy_solution failed");

        STATS_increment_added_cuts();
    }

CLEANUP:
    if (x) free(x);
    return rval;
}

static void dump_row(FILE *file, const struct Row *row, int k)
{
    fprintf(file, "double pi%d[] = { ", k);
    for(int i = 0; i < row->nz; i++)
        fprintf(file, "%.10lf, ", row->pi[i]);
    fprintf(file, "};\n");

    fprintf(file, "int indices%d[] = { ", k);
    for(int i = 0; i < row->nz; i++)
        fprintf(file, "%d, ", row->indices[i]);
    fprintf(file, "};\n");

    fprintf(file, "struct Row row%d = {", k);
    fprintf(file, ".nz = %d, ", row->nz);
    fprintf(file, ".head = %d, ", row->head);
    fprintf(file, ".pi_zero = %.10lf, ", row->pi_zero);
    fprintf(file, ".pi = pi%d, ", k);
    fprintf(file, ".indices = indices%d };\n", k);
}

static void dump_cut(const struct Row *cut, long count)
{
    char filename[100];
    sprintf(filename, "cut-%03ld.c", count);

    log_info("Dumping generated cut to file %s\n", filename);

    FILE *file = fopen(filename, "w");
    dump_row(file, cut, 0);
    fclose(file);
}

static void dump_tableau(const struct Tableau *tableau, long count)
{
    char filename[100];
    sprintf(filename, "tableau-%03ld.c", count);

    log_info("Dumping tableau to file %s\n", filename);

    FILE *file = fopen(filename, "w");

    for(int i = 0; i < tableau->nrows; i++)
        dump_row(file, tableau->rows[i], i);

    fprintf(file, "struct Row* rows[] = {");
    for(int i = 0; i < tableau->nrows; i++)
        fprintf(file, "&row%d, ", i);
    fprintf(file, "};\n");

    fprintf(file, "char column_types[] = {");
    for(int i = 0; i < tableau->ncols; i++)
        fprintf(file, "%d, ", tableau->column_types[i]);
    fprintf(file, "};\n");

    fprintf(file, "struct Tableau tableau = {");
    fprintf(file, ".ncols = %d, ", tableau->ncols);
    fprintf(file, ".nrows = %d, ", tableau->nrows);
    fprintf(file, ".rows = rows, ");
    fprintf(file, ".column_types = column_types };\n");

    fclose(file);
}

static void dump_integral_solution(const struct CG *cg, long count)
{
    char filename[100];
    sprintf(filename, "solution-%03ld.c", count);

    log_info("Dumping integral solution to file %s\n", filename);

    FILE *file = fopen(filename, "w");
    double *x = cg->integral_solution;

    fprintf(file, "x[%d] = {", cg->ncols);
    for(int i = 0; i < cg->ncols; i++)
        fprintf(file, "%.10lf, ", x[i]);
    fprintf(file, "};\n");

    fclose(file);
}

#ifndef TEST_SOURCE

double CG_replace_x(const struct Row *row, const double *x)
{
    double lhs = 0;
    double *terms = (double*) malloc(row->nz * sizeof(double));

    for (int i = 0; i < row->nz; i++)
    {
        double pii = row->pi[i];
        int idx = row->indices[i];
        terms[i] = pii * x[idx];
    }

    DOUBLE_sum(terms, row->nz, &lhs);

    free(terms);
    return lhs;
}

/*
 * Looks for a certain ray at an array of rays. Matches rays that are
 * parallel to the given one. Returns whether a matching ray was found,
 * the index of the matching ray and the scaling factor.
 */
int CG_find_ray(const struct RayList *rays,
                const double *r,
                int *found,
                double *scale,
                int *index)
{
    *found = 0;

    for (int i = 0; i < rays->nrays; i++)
    {
        double *q = LFREE_get_ray(rays, i);

        int match;
        check_rays_parallel(rays->dim, r, q, &match, scale);

        if (match)
        {
            *index = i;
            *found = 1;
            goto CLEANUP;
        }
    }

CLEANUP:
    return 0;
}

int CG_extract_model(const struct Tableau *tableau,
                     struct TableauModelMap *map,
                     struct MultiRowModel *model)
{
    int rval = 0;

    int nrows = tableau->nrows;
    struct Row **rows = tableau->rows;

    map->nvars = 0;

    int *i = 0;
    int *idx = 0;

    i = (int *) malloc(nrows * sizeof(int));
    idx = (int *) malloc(nrows * sizeof(int));
    abort_if(!i, "could not allocate i");
    abort_if(!idx, "could not allocate idx");

    for (int j = 0; j < nrows; j++)
        i[j] = 0;

    while (1)
    {
        double *r = LFREE_get_ray(&model->rays, model->rays.nrays);

        int idx_min = INT_MAX;

        for (int j = 0; j < nrows; j++)
        {
            r[j] = 0.0;

            if (i[j] < rows[j]->nz)
                idx[j] = rows[j]->indices[i[j]];
            else
                idx[j] = INT_MAX;

            idx_min = min(idx_min, idx[j]);
        }

        if (idx_min == INT_MAX)
            break;

        for (int j = 0; j < nrows; j++)
        {
            if (idx[j] > idx_min) continue;
            r[j] = -rows[j]->pi[i[j]];
            i[j]++;
        }

        for (int j = 0; j < nrows; j++)
            if (idx_min == rows[j]->head)
                goto NEXT_RAY;

        int found;
        double scale;
        int ray_index;

        log_verbose("  extracted ray (%d):\n", idx_min);
        for (int j = 0; j < nrows; j++)
                log_verbose("    r[%d] = %.12lf\n", j, r[j]);

        rval = CG_find_ray(&model->rays, r, &found, &scale, &ray_index);
        abort_if(rval, "CG_find_ray failed");

        found = 0;
        if (!found)
        {
            log_verbose("  ray is new\n");
            scale = 1.0;
            ray_index = model->rays.nrays++;
        }
        else
        {
            log_verbose("  ray equals:\n");
            double *q = LFREE_get_ray(&model->rays, ray_index);
            for (int j = 0; j < nrows; j++)
                log_verbose("    r[%d] = %.12lf\n", j, q[j]);
            log_verbose("    scale = %.12lf\n", scale);
        }

        map->ray_scale[map->nvars] = scale;
        map->indices[map->nvars] = idx_min;
        map->variable_to_ray[map->nvars] = ray_index;
        map->nvars++;

    NEXT_RAY:;
    }

//    for (int j = 0; j < model->rays.nrays; j++)
//    {
//        double *r = LFREE_get_ray(&model->rays, j);
//        double max_scale = 0.0;
//
//        for (int k = 0; k < map->nvars; k++)
//        {
//            if (map->variable_to_ray[k] != j) continue;
//            if (map->ray_scale[k] < max_scale) continue;
//            max_scale = map->ray_scale[k];
//        }
//
//        abort_if(max_scale == 0.0, "max_scale is zero");
//
//        for (int k = 0; k < map->nvars; k++)
//            if (map->variable_to_ray[k] == j)
//                map->ray_scale[k] /= max_scale;
//
//        for (int k = 0; k < nrows; k++)
//            r[k] *= max_scale;
//    }

CLEANUP:
    if (idx) free(idx);
    if (i) free(i);
    return rval;
}

int CG_init(struct LP *lp, char *column_types, struct CG *cg)
{
    int rval = 0;

    cg->lp = lp;
    cg->column_types = column_types;

    cg->tableau_rows = 0;
    cg->cstat = 0;
    cg->rstat = 0;
    cg->ub = 0;
    cg->lb = 0;
    cg->last_obj_value = NAN;

    cg->integral_solution = 0;
    cg->basic_solution = 0;
    cg->current_solution = 0;

    int nrows = LP_get_num_rows(lp);
    int ncols = LP_get_num_cols(lp);

    cg->nrows = nrows;
    cg->ncols = ncols;

    cg->cstat = (int *) malloc(ncols * sizeof(int));
    cg->rstat = (int *) malloc(nrows * sizeof(int));
    cg->ub = (double *) malloc(ncols * sizeof(double));
    cg->lb = (double *) malloc(ncols * sizeof(double));

    abort_if(!cg->cstat, "could not allocate cg->cstat");
    abort_if(!cg->rstat, "could not allocate cg->rstat");
    abort_if(!cg->ub, "could not allocate cg->ub");
    abort_if(!cg->lb, "could not allocate cg->lb");

    cg->tableau_rows = (struct Row **) malloc(nrows * sizeof(struct Row *));
    abort_if(!cg->tableau_rows, "could not allocate cg->tableau_rows");

    rval = LP_get_tableau(lp, cg->tableau_rows, cg->cstat, cg->rstat, cg->ub,
            cg->lb);
    abort_if(rval, "LP_get_tableau failed");

    cg_initial_time = get_user_time();

CLEANUP:
    return rval;
}

void CG_free(struct CG *cg)
{
    if (!cg) return;
    if (cg->ub) free(cg->ub);
    if (cg->lb) free(cg->lb);
    if (cg->cstat) free(cg->cstat);
    if (cg->rstat) free(cg->rstat);

    if (cg->integral_solution) free(cg->integral_solution);
    if (cg->basic_solution) free(cg->basic_solution);
    if (cg->current_solution) free(cg->current_solution);

    if (cg->tableau_rows)
    {
        for (int i = 0; i < cg->nrows; i++)
        {
            LP_free_row(cg->tableau_rows[i]);
            free(cg->tableau_rows[i]);
        }
        free(cg->tableau_rows);
    }

    free(cg);
}

int CG_add_single_row_cuts(struct CG *cg, SingleRowGeneratorCallback generate)
{
    int rval = 0;

    for (int i = 0; i < cg->nrows; i++)
    {
        struct Row *row = cg->tableau_rows[i];

        if (frac(row->pi_zero) < EPSILON)
            continue;

        if (cg->column_types[row->head] != MILP_INTEGER)
            continue;

        log_verbose("Generating cut %d...\n", i);

        struct Row cut;

        rval = generate(row, cg->column_types, &cut);
        abort_if(rval, "generate failed");

        int ignored;
        rval = add_cut(cg, &cut, &ignored);
        abort_if(rval, "add_cut failed");

        LP_free_row(&cut);
    }

CLEANUP:
    return rval;
}

int CG_estimate_multirow_cut_count(struct CG *cg,
                                   int nrows,
                                   int *row_selected,
                                   int *row_affinity,
                                   long *total_count,
                                   double score_cutoff)
{
    int rval = 0;
    int *row_indices = 0;
    int finished = 0;

    *total_count = 0;

    row_indices = (int *) malloc(nrows * sizeof(int));
    abort_if(!row_indices, "could not allocate row_indices");

    for (int i = 0; i < nrows; i++)
        row_indices[i] = nrows - i - 1;

    for (int i = 0; i < cg->nrows * cg->nrows; i++)
        row_affinity[i] = -1;

    do
    {
        int inc_index = 0;
        int is_rhs_integer = 1;
        int valid_combination = 1;

        for (int i = 0; i < nrows; i++)
        {
            struct Row *r = cg->tableau_rows[row_indices[i]];

            if (!row_selected[row_indices[i]])
            {
                valid_combination = 0;
                inc_index = i;
            }

            double df = fabs(frac(r->pi_zero) - 0.5);
            if (df < INTEGRALITY_THRESHOLD) is_rhs_integer = 0;
        }

        if (is_rhs_integer) valid_combination = 0;

        for (int i = 0; valid_combination && i < nrows; i++)
        {
            for (int j = i + 1; valid_combination && j < nrows; j++)
            {
                int i1 = row_indices[i];
                int i2 = row_indices[j];
                if (i2 < i1) swap(i1, i2, int);

                int k = cg->nrows * i1 + i2;

                if (row_affinity[k] < 0)
                {
                    struct Row *row1 = cg->tableau_rows[i1];
                    struct Row *row2 = cg->tableau_rows[i2];
                    double score;

                    rval = evaluate_row_pair(row1, row2, &score);
                    abort_if(rval, "evaluate_row_pair failed");

                    row_affinity[k] = 1;
                    if (score < score_cutoff) row_affinity[k] = 0;
                }

                if (!row_affinity[k])
                {
                    valid_combination = 0;
                    goto NEXT_COMBINATION;
                }
            }
        }

        if (valid_combination)
        {
            (*total_count)++;
        }

    NEXT_COMBINATION:
        next_combination(cg->nrows, nrows, row_indices, inc_index, &finished);
    } while (!finished);

CLEANUP:
    if (row_indices) free(row_indices);
    return rval;
}

int CG_cut_dynamism(struct Row *cut, double *dynamism)
{
    double largest = -INFINITY;
    double smallest = INFINITY;

    for (int i = 0; i < cut->nz; i++)
    {
        smallest = fmin(fabs(cut->pi[i]), smallest);
        largest = fmax(fabs(cut->pi[i]), largest);
    }

    *dynamism = largest / smallest;

    return 0;
}

int CG_add_multirow_cuts(struct CG *cg,
                         int nrows,
                         MultiRowGeneratorCallback generate)
{
    int rval = 0;
    int *row_indices = 0;
    struct Row **rows = 0;
    int *row_affinity = 0;
    int *row_selected = 0;

    int finished;
    long count = 0;
    long total_count = 0;

    row_indices = (int *) malloc(nrows * sizeof(int));
    rows = (struct Row **) malloc(nrows * sizeof(struct Row *));
    row_affinity = (int *) malloc((cg->nrows * cg->nrows) * sizeof(int));
    row_selected = (int *) malloc(cg->nrows * sizeof(int));

    abort_if(!row_indices, "could not allocate row_indices");
    abort_if(!rows, "could not allocate rows");
    abort_if(!row_affinity, "could not allocate row_affinity");
    abort_if(!row_selected, "could not allocate row_selected");

    rval = select_rows(cg, row_selected);
    abort_if(rval, "select_rows failed");

    log_info("    Finding combinations...\n");
    for (double cutoff = 0.05; cutoff <= 1.0; cutoff += 0.05)
    {
        double cg_current_time = get_user_time() - cg_initial_time;
        if (cg_current_time > CG_TIMEOUT) break;

        rval = CG_estimate_multirow_cut_count(cg, nrows, row_selected,
                row_affinity, &total_count, cutoff);
        abort_if(rval, "estimate_two_row_cuts_count failed");

        log_verbose("    %8d combinations [%.2lf]\n", total_count, cutoff);

        if (total_count < MAX_SELECTED_COMBINATIONS)
        {
            log_info("    %8d combinations [%.2lf]\n", total_count, cutoff);
            break;
        }
    }

    for (int i = 0; i < nrows; i++)
        row_indices[i] = nrows - i - 1;

    total_count = min(total_count, MAX_SELECTED_COMBINATIONS);

    progress_set_total(total_count);
    progress_reset();

    do
    {
        double cg_current_time = get_user_time() - cg_initial_time;
        if (cg_current_time > CG_TIMEOUT) break;

        int inc_index = 0;
        int is_rhs_integer = 1;
        int valid_combination = 1;

        for (int i = 0; i < nrows; i++)
        {
            rows[i] = cg->tableau_rows[row_indices[i]];

            if (!row_selected[row_indices[i]])
            {
                valid_combination = 0;
                inc_index = i;
            }

            double df = fabs(frac(rows[i]->pi_zero) - 0.5);
            if (df < INTEGRALITY_THRESHOLD) is_rhs_integer = 0;
        }

        if (is_rhs_integer) valid_combination = 0;

        for (int i = 0; valid_combination && i < nrows; i++)
        {

            for (int j = i + 1; valid_combination && j < nrows; j++)
            {
                int i1 = row_indices[i];
                int i2 = row_indices[j];
                if (i2 < i1) swap(i1, i2, int);

                int k = cg->nrows * i1 + i2;

                abort_if(row_affinity[k] < 0, "row_affinity not computed");
                if (!row_affinity[k]) valid_combination = 0;

                if_verbose_level
                {
                    struct Row *row1 = cg->tableau_rows[i1];
                    struct Row *row2 = cg->tableau_rows[i2];
                    double score;

                    rval = evaluate_row_pair(row1, row2, &score);
                    abort_if(rval, "evaluate_row_pair failed");

                    log_verbose("%4d %4d %.2lf\n", i1, i2, score);
                }
            }
        }

        if (valid_combination)
        {
            count++;
            if (count > MAX_SELECTED_COMBINATIONS)
            {
                log_info(
                        "    maximum number of combinations reached. stopping.\n");
                break;
            }

            if_debug_level if (ONLY_CUT > 0 && count != ONLY_CUT)
                    goto NEXT_COMBINATION;

            if_verbose_level
            {
                time_printf("Generating cut %d from [ ", count);
                for (int i = 0; i < nrows; i++)
                    printf("%d ", row_indices[i]);
                printf("]...\n");
            }

            if (LOG_LEVEL == LOG_LEVEL_INFO)
            {
                progress_print();
                progress_increment();
            }

            struct Tableau tableau =
            {
                    .ncols = cg->ncols,
                    .nrows = nrows,
                    .rows = rows,
                    .column_types = cg->column_types
            };

            if_verbose_level dump_tableau(&tableau, count);

            int max_nz = CG_total_nz(&tableau);

            struct Row cut;
            rval = LP_init_row(&cut, max_nz);
            abort_if(rval, "LP_init_row failed");

            double initial_time = get_user_time();

            rval = generate(&tableau, &cut);
            if (rval == ERR_NO_CUT)
            {
                rval = 0;
                log_verbose("combination does not yield cut\n");
                LP_free_row(&cut);
                goto NEXT_COMBINATION;
            }
            else if(rval) {
                dump_tableau(&tableau, count);
                abort_iff(1, "generate failed (cut %d)", count);
            }

            if_verbose_level dump_cut(&cut, count);

            double elapsed_time = get_user_time() - initial_time;
            log_verbose("    generate: %.2lf ms\n", elapsed_time * 1000);

            double dynamism;
            rval = CG_cut_dynamism(&cut, &dynamism);
            abort_if(rval, "CG_cut_dynamism failed");

            if (dynamism > MAX_CUT_DYNAMISM)
            {
                log_verbose("Discarding cut (dynamism=%.2lf)\n", dynamism);
                LP_free_row(&cut);
                goto NEXT_COMBINATION;
            }

            int ignored;
            rval = add_cut(cg, &cut, &ignored);
            if(rval)
            {
                dump_cut(&cut, count);
                dump_tableau(&tableau, count);
                dump_integral_solution(cg, count);
                abort_iff(1, "add_cut failed (cut=%d)", count);
            }

            if_debug_level if (!ignored)
            {
                SHOULD_DUMP_CUTS = 1;
                generate(&tableau, &cut);
                SHOULD_DUMP_CUTS = 0;
            }

            LP_free_row(&cut);
        }

    NEXT_COMBINATION:
        next_combination(cg->nrows, nrows, row_indices, inc_index, &finished);
    } while (!finished);

CLEANUP:
    if (row_selected) free(row_selected);
    if (row_affinity) free(row_affinity);
    if (row_indices) free(row_indices);
    if (rows) free(rows);
    return rval;
}

int CG_set_integral_solution(struct CG *cg, double *valid_solution)
{
    return copy_solution(cg, valid_solution, &cg->integral_solution);
}

int CG_set_basic_solution(struct CG *cg, double *basic_solution)
{
    int rval = 0;

    rval = copy_solution(cg, basic_solution, &cg->basic_solution);
    abort_if(rval, "copy_solution failed");

    rval = copy_solution(cg, basic_solution, &cg->current_solution);
    abort_if(rval, "copy_solution failed");

CLEANUP:
    return rval;
}

int CG_boost_variable(int var,
                      double factor,
                      int nrows,
                      double *rays,
                      int *variable_to_ray,
                      double *ray_scale,
                      int *indices,
                      int nz)
{
    int rval = 0;

    int selected_ray = -1;
    for (int j = 0; j < nz; j++)
    {
        if (indices[j] == var)
        {
            selected_ray = variable_to_ray[j];
            break;
        }
    }

    if (selected_ray < 0) goto CLEANUP;

    double *r = &rays[nrows * selected_ray];
    for (int k = 0; k < nrows; k++)
        r[k] *= factor;

    for (int j = 0; j < nz; j++)
    {
        if (variable_to_ray[j] == selected_ray)
            ray_scale[j] /= factor;
    }

CLEANUP:
    return rval;
}

int CG_init_map(struct TableauModelMap *map, int max_nrays, int nrows)
{
    int rval = 0;

    map->variable_to_ray = (int *) malloc(max_nrays * sizeof(int));
    map->indices = (int *) malloc(max_nrays * sizeof(int));
    map->ray_scale = (double *) malloc(max_nrays * sizeof(double));
    abort_if(!map->variable_to_ray, "could not allocate variable_to_ray");
    abort_if(!map->indices, "could not allocate indices");
    abort_if(!map->ray_scale, "could not allocate ray_scale");

CLEANUP:
    return rval;
}

void CG_free_map(struct TableauModelMap *map)
{
    if(!map) return;
    free(map->variable_to_ray);
    free(map->indices);
    free(map->ray_scale);
}

int CG_extract_f_from_tableau(const struct Tableau *tableau, double *f)
{
    for (int i = 0; i < tableau->nrows; i++)
    {
        f[i] = frac(tableau->rows[i]->pi_zero);
        if (DOUBLE_eq(f[i], 1.0)) f[i] = 0.0;
    }

    return 0;
}

int CG_init_model(struct MultiRowModel *model, int nrows, int rays_capacity)
{
    int rval = 0;
    model->nrows = nrows;

    rval = LFREE_init_ray_list(&model->rays, nrows, rays_capacity);
    abort_if(rval, "LFREE_init_ray_list failed");

    model->f = (double*) malloc(nrows * sizeof(double));
    abort_if(!model->f, "could not allocate f");

CLEANUP:
    return rval;
}

void CG_free_model(struct MultiRowModel *model)
{
    if(!model) return;
    free(model->f);
    LFREE_free_ray_list(&model->rays);
}

int CG_total_nz(const struct Tableau *tableau)
{
    int total_nz = 0;
    for(int i = 0; i < tableau->nrows; i++)
        total_nz += tableau->rows[i]->nz;
    return total_nz;
}

int CG_print_model(const struct MultiRowModel *model)
{
    int rval = 0;
    int nrows = model->nrows;

    time_printf("       f = [");
    for (int i = 0; i < nrows; i++)
        printf("%20.12lf ", model->f[i]);
    printf("]\n");

    for (int i = 0; i < model->rays.nrays; i++)
    {
        double *ray = LFREE_get_ray(&model->rays, i);
        double norm = 0;

        time_printf("ray[%3d] = [", i);
        for (int j = 0; j < nrows; j++)
        {
            printf("%20.12lf ", ray[j]);
            norm += fabs(ray[j]);
        }
        printf("] norm=%20.12lf \n", norm);
    }

CLEANUP:
    return rval;
}

#endif // TEST_SOURCE
