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

#include <infinity/infinity.h>
#include <infinity/infinity-2d.h>
#include <infinity/greedy-nd.h>

struct SortPair
{
    int index;
    void *data;
};

/**
 * Compares two rays according to their angles.
 *
 * The input are two SortPairs, where the data points to a double[2]. If the
 * two rays are pointing to exactly opposite direction, returns 1.
 *
 * @param p1 a pointer to a SortPair containing the first ray.
 * @param p2 a pointer to a SortPair containing the second ray.
 * @return -1, 1 or 0, if the second ray is in clockwise, counter-clockwise
 *         or aligned, respectively, to the first ray.
 */
static int _qsort_cmp_rays_angle(const void *p1, const void *p2)
{
    double *r1 = (double *) (((struct SortPair *) p1)->data);
    double *r2 = (double *) (((struct SortPair *) p2)->data);
    return sign(atan2(r1[0], r1[1]) - atan2(r2[0], r2[1]));
}

/**
 * Sorts a list of rays according to their angle.
 *
 * @param rays the rays to be sorted
 * @param nrays the number of rays in the list
 * @param beta a list of doubles, to be sorted together with the rays
 *
 * @return zero if successful, non-zero otherwise
 */
static int sort_rays_angle(double *rays, int nrays, double *beta)
{
    int rval = 0;

    double *rays_copy = 0;
    double *beta_copy = 0;
    struct SortPair *pairs = 0;

    pairs = (struct SortPair *) malloc(nrays * sizeof(struct SortPair));
    rays_copy = (double *) malloc(2 * nrays * sizeof(double));
    beta_copy = (double *) malloc(nrays * sizeof(double));
    abort_if(!pairs, "could not allocate pairs");
    abort_if(!rays_copy, "could not allocate rays_copy");
    abort_if(!beta_copy, "could not allocate beta_copy");

    memcpy(rays_copy, rays, 2 * nrays * sizeof(double));
    memcpy(beta_copy, beta, nrays * sizeof(double));

    for (int i = 0; i < nrays; i++)
    {
        pairs[i].index = i;
        pairs[i].data = &rays[2 * i];
    }

    qsort(pairs, (size_t) nrays, sizeof(struct SortPair),
            _qsort_cmp_rays_angle);

    for (int i = 0; i < nrays; i++)
    {
        beta[i] = beta_copy[pairs[i].index];
        memcpy(&rays[2 * i], &rays_copy[2 * pairs[i].index],
                2 * sizeof(double));
    }

CLEANUP:
    if (pairs) free(pairs);
    if (rays_copy) free(rays_copy);
    if (beta_copy) free(beta_copy);
    return rval;
}

static void print_row(const struct Row *row)
{
    time_printf("Row:\n");
    for (int i = 0; i < row->nz; i++)
        time_printf("    %.4lfx%d\n", row->pi[i], row->indices[i]);
    time_printf("    <= %.4lf [%d]\n", row->pi_zero, row->head);
}

static void print_rays(const struct Tableau *tableau,
                       const struct RayMap *map,
                       const double *f,
                       const double *rays,
                       int nrays)
{
    int nrows = tableau->nrows;

    time_printf("Ray map:\n");
    for (int i = 0; i < map->nvars; i++)
        time_printf("    %4d: %4d x%-4d (scale=%.4lf)\n", i,
                map->variable_to_ray[i], map->indices[i], map->ray_scale[i]);

    time_printf("Origin:\n");
    for (int i = 0; i < nrows; i++)
        time_printf("    %20.12lf\n", f[i]);

    time_printf("Rays:\n");
    for (int i = 0; i < nrays; i++)
    {
        time_printf("    ");

        for (int j = 0; j < nrows; j++)
            printf("%20.12lf ", rays[i * nrows + j]);

        printf("angle=%.4lf ", atan2(rays[i * nrows], rays[i * nrows + 1]));
        printf("norm=%.4lf ",
                fabs(rays[i * nrows]) + fabs(rays[i * nrows + 1]));

        printf("[ ");
        for (int j = 0; j < map->nvars; j++)
            if (map->variable_to_ray[j] == i)
                printf("%d ", map->indices[j]);
        printf("]\n");
    }
}

static void print_cut(const struct Row *cut)
{
    time_printf("Generated cut:\n");
    for (int i = 0; i < cut->nz; i++)
        time_printf("    %.4lfx%d\n", cut->pi[i], cut->indices[i]);
    time_printf("     <= %.4lf\n", cut->pi_zero);
}

static int create_cut(const struct Tableau *tableau,
                      const struct RayMap *map,
                      const double *f,
                      const int lfree_nrays,
                      const double *lfree_rays,
                      const double *beta,
                      struct Row *cut)
{
    int rval = 0;
    int nvars = map->nvars;
    int nrows = tableau->nrows;

    cut->nz = nvars;
    cut->pi = (double *) malloc(nvars * sizeof(double));
    cut->indices = (int *) malloc(nvars * sizeof(int));
    abort_if(!cut->pi, "could not allocate cut->pi");
    abort_if(!cut->indices, "could not allocate cut->indices");

    struct LP lp;

    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    rval = GREEDY_create_psi_lp(nrows, lfree_nrays, f, lfree_rays, beta, &lp);
    abort_if(rval, "create_psi_lp failed");

    for (int i = 0; i < nvars; i++)
    {
        double value;
        const double *q = &map->rays[map->variable_to_ray[i] * nrows];

        if (ENABLE_LIFTING &&
                tableau->column_types[map->indices[i]] == MILP_INTEGER)
        {
            rval = GREEDY_ND_pi(nrows, lfree_nrays, f, lfree_rays, beta, q,
                    map->ray_scale[i], &lp, &value);
            abort_if(rval, "GREEDY_ND_pi failed");
        }
        else
        {
            rval = GREEDY_ND_psi(nrows, lfree_nrays, f, lfree_rays, beta, q,
                    map->ray_scale[i], &lp, &value);
            abort_if(rval, "GREEDY_ND_psi failed");
        }

        log_verbose("   psi[%4d] = %20.12lf %d\n", map->indices[i], value);

        value *= 1.0001;
        value = DOUBLE_max(value, 0.0001);

        cut->indices[i] = map->indices[i];
        cut->pi[i] = -value;
    }

    cut->pi_zero = -1.0;

CLEANUP:
    LP_free(&lp);
    return rval;
}

static int select_rays(const struct RayMap map,
                       const struct Tableau *tableau,
                       struct MultiRowModel *model)
{
    int rval = 0;
    int nrows = tableau->nrows;

    for (double norm_cutoff = 0.00; norm_cutoff <= 5.0; norm_cutoff += 0.1)
    {
        model->nrays = 0;

        for (int i = 0; i < map.nrays; i++)
        {
            int keep = 1;
            double *r = &map.rays[tableau->nrows * i];

            for (int j = 0; j < (model->nrays); j++)
            {
                double *q = &map.rays[tableau->nrows * j];
                double norm = 0;

                for (int k = 0; k < nrows; k++)
                    norm += fabs(r[k] - q[k]);

                if (norm <= norm_cutoff)
                {
                    keep = 0;
                    break;
                }
            }

            if (keep)
            {
                memcpy(&model->rays[nrows * (model->nrays)], r,
                        nrows * sizeof(double));
                model->nrays++;
            }
        }

        log_debug("  norm_cutoff=%8.2lf nrays=%8d\n", norm_cutoff,
                model->nrays);

        if (model->nrays < MAX_N_RAYS) break;
    }

CLEANUP:
    return rval;
}

static int
append_extra_rays(const struct Tableau *tableau, double *rays, int *nrays)
{
    int rval = 0;
    int nrows = tableau->nrows;
    int n_extra_rays = 0;
    double extra_rays[100];

    abort_if(nrows > 3, "not implemented");

    if (nrows == 2)
    {
        extra_rays[0] = 0.0001;
        extra_rays[1] = 0.0000;

        extra_rays[2] = -0.0001;
        extra_rays[3] = 0.0001;

        extra_rays[4] = -0.0001;
        extra_rays[5] = -0.0001;

        n_extra_rays = 3;
    }

    if (nrows == 3)
    {
        extra_rays[0] = 0.0000;
        extra_rays[1] = 0.0000;
        extra_rays[2] = 0.0001;

        extra_rays[3] = 0.0001;
        extra_rays[4] = 0.0000;
        extra_rays[5] = -0.0001;

        extra_rays[6] = -0.0001;
        extra_rays[7] = 0.0000;
        extra_rays[8] = -0.0001;

        extra_rays[9] = -0.0001;
        extra_rays[10] = -0.0001;
        extra_rays[11] = -0.0001;

        n_extra_rays = 4;
    }

    for (int i = 0; i < n_extra_rays; i++)
    {
        double *r = &extra_rays[nrows * i];

        double scale;
        int found, index;
        rval = CG_find_ray(nrows, rays, *nrays, r, &found, &scale, &index);
        abort_if(rval, "CG_find_ray failed");

        if (!found)
        {
            memcpy(&rays[nrows * (*nrays)], r, nrows *
                    sizeof(double));
            (*nrays)++;
        }
    }

CLEANUP:
    return rval;
}

static int write_sage_file(int nrows,
                           int nrays,
                           const double *f,
                           const double *rays,
                           const double *beta,
                           const char *filename)
{
    int rval = 0;

    FILE *fsage = fopen(filename, "w");
    abort_iff(!fsage, "could not open %s", filename);

    fprintf(fsage, "f=vector([");
    for (int i = 0; i < nrows; i++)
        fprintf(fsage, "%.20lf,", f[i]);
    fprintf(fsage, "])\n");

    fprintf(fsage, "R=matrix([\n");
    for (int i = 0; i < nrays; i++)
    {
        fprintf(fsage, "    [");
        for (int j = 0; j < nrows; j++)
        {
            fprintf(fsage, "%.20lf,", rays[i * nrows + j]);
        }
        fprintf(fsage, "],\n");
    }
    fprintf(fsage, "])\n");

    fprintf(fsage, "pi=vector([\n");
    for (int k = 0; k < nrays; k++)
        fprintf(fsage, "    %.12lf,\n", 1 / beta[k]);
    fprintf(fsage, "])\n");

CLEANUP:
    fclose(fsage);
    return rval;
}

static int dump_cut(const struct Tableau *tableau,
                    const double *rays,
                    int nrays,
                    const double *f,
                    const double *beta)
{
    int rval = 0;

    char filename[100];
    sprintf(filename, "cut-%03d.sage", DUMP_CUT_N++);

    time_printf("Writing %s...\n", filename);
    rval = write_sage_file(tableau->nrows, nrays, f, rays, beta, filename);
    abort_if(rval, "write_sage_file failed");

CLEANUP:
    return rval;
}

#ifndef TEST_SOURCE

int INFINITY_generate_cut(struct Tableau *tableau, struct Row *cut)
{
    int rval = 0;
    int max_nrays = 0;
    int nrows = tableau->nrows;
    double *beta = 0;

    for (int i = 0; i < tableau->nrows; i++)
        max_nrays += tableau->rows[i]->nz;

    struct MultiRowModel model;
    rval = CG_init_model(&model, nrows, max_nrays + 100);
    abort_if(rval, "CG_init_model failed");

    rval = CG_extract_f_from_tableau(tableau, model.f);
    abort_if(rval, "CG_extract_f_from_tableau failed");

    struct RayMap map;
    rval = CG_init_ray_map(&map, max_nrays, nrows);
    abort_if(rval, "CG_init_ray_map failed");

    rval = CG_extract_rays_from_tableau(tableau, &map);
    abort_if(rval, "CG_extract_rays_from_rows failed");

    rval = select_rays(map, tableau, &model);
    abort_if(rval, "select_rays failed");

    if (ENABLE_LIFTING)
    {
        rval = append_extra_rays(tableau, model.rays, &model.nrays);
        abort_if(rval, "append_extra_rays failed");
    }

    beta = (double *) malloc(model.nrays * sizeof(double));
    abort_if(!beta, "could not allocate beta");

    if (model.nrays < 3)
    {
        rval = ERR_NO_CUT;
        cut->pi = 0;
        cut->indices = 0;
        goto CLEANUP;
    }

    log_debug("Selected %d rays\n", nrays);
    if_verbose_level print_rays(tableau, &map, model.f, map.rays, map.nrays);

    log_verbose("Computing lattice-free set...\n");
    if (nrows == 2)
    {
        rval = sort_rays_angle(model.rays, model.nrays, beta);
        abort_if(rval, "sort_rays_angle failed");

        rval = INFINITY_2D_generate_cut(&model, beta);
        if (rval)
        {
            rval = ERR_NO_CUT;
            goto CLEANUP;
        }
    }
    else
    {
        rval = GREEDY_ND_generate_cut(nrows, model.nrays, model.f, model.rays,
                beta);
        if (rval)
        {
            rval = ERR_NO_CUT;
            goto CLEANUP;
        }
    }

    if (SHOULD_DUMP_CUTS)
    {
        rval = dump_cut(tableau, model.rays, model.nrays, model.f, beta);
        abort_if(rval, "dump_cut failed");
    }

    rval = create_cut(tableau, &map, model.f, model.nrays, model.rays, beta,
            cut);
    abort_if(rval, "create_cut failed");

    if_verbose_level print_cut(cut);

CLEANUP:
    if (beta) free(beta);
    CG_free_model(&model);
    CG_free_ray_map(&map);
    return rval;
}

#endif // TEST_SOURCE
