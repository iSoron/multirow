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

#include <infinity/greedy.h>
#include <infinity/greedy-bsearch.h>
#include <infinity/greedy-2d.h>
#include <infinity/greedy-nd.h>

struct SortPair
{
    int index;
    void *data;
};

static int _qsort_cmp_rays_angle(const void *p1,
                                 const void *p2)
{
    double *r1 = (double *) (((struct SortPair *) p1)->data);
    double *r2 = (double *) (((struct SortPair *) p2)->data);
    return sign(atan2(r1[0], r1[1]) - atan2(r2[0], r2[1]));
}

static int sort_rays_angle(double *rays,
                           int nrays,
                           double *beta)
{
    int rval = 0;

    int *map = 0;
    double *rays_copy = 0;
    double *beta_copy = 0;
    struct SortPair *pairs = 0;

    pairs = (struct SortPair *) malloc(nrays * sizeof(struct SortPair));
    rays_copy = (double *) malloc(2 * nrays * sizeof(double));
    beta_copy = (double*) malloc(nrays * sizeof(double));
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

    qsort(pairs, nrays, sizeof(struct SortPair), _qsort_cmp_rays_angle);

    for (int i = 0; i < nrays; i++)
    {
        beta[i] = beta_copy[pairs[i].index];
        memcpy(&rays[2 * i], &rays_copy[2 * pairs[i].index], 2 *
                sizeof(double));
    }

CLEANUP:
    if (pairs) free(pairs);
    if (rays_copy) free(rays_copy);
    if (beta_copy) free(beta_copy);
    return rval;
}


static int scale_rays(double *rays,
                      int nrays,
                      double *scale)
{
    int rval = 0;

    for (int i = 0; i < nrays; i++)
    {
        rays[2*i] *= scale[i];
        rays[2*i+1] *= scale[i];
    }

CLEANUP:
    return rval;
}


static void print_row(const struct Row *row)
{
    time_printf("Row:\n");
    for (int i = 0; i < row->nz; i++)
        time_printf("    %.4lfx%d\n", row->pi[i], row->indices[i]);
    time_printf("    <= %.4lf [%d]\n", row->pi_zero, row->head);
}

static void print_rays(const int *map,
                       const int *indices,
                       const double *f,
                       const double *rays,
                       const double *scale,
                       int nrays,
                       int nz,
                       int nrows)
{
    time_printf("Mapping:\n");
    for (int i = 0; i < nz; i++)
        time_printf("    %4d: %4d x%-4d (scale=%.4lf)\n", i, map[i], indices[i],
                    scale[i]);

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
        printf("norm=%.4lf ", fabs(rays[i * nrows]) + fabs(rays[i * nrows + 1]));

        printf("[ ");
        for (int j = 0; j < nz; j++)
            if (map[j] == i) printf("%d ", indices[j]);
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

int GREEDY_write_sage_file(int nrows,
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

static int create_cut(const int nrows,
                      const char *column_types,
                      const int nrays,
                      const double *f,
                      const int lfree_nrays,
                      const double *lfree_rays,
                      const double *beta,
                      const double *extracted_rays,
                      const int nvars,
                      const int *variable_to_ray,
                      const int *indices,
                      const double *scale,
                      struct Row *cut)
{
    int rval = 0;

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
        const double *q = &extracted_rays[variable_to_ray[i] * nrows];

        if(ENABLE_LIFTING && column_types[indices[i]] == MILP_INTEGER)
        {
            rval = GREEDY_ND_pi(nrows, lfree_nrays, f, lfree_rays, beta, q,
                    scale[i], &lp, &value);
            abort_if(rval, "GREEDY_ND_pi failed");
        }
        else
        {
            rval = GREEDY_ND_psi(nrows, lfree_nrays, f, lfree_rays, beta, q,
                    scale[i], &lp, &value);
            abort_if(rval, "GREEDY_ND_psi failed");
        }

        log_verbose("   psi[%4d] = %20.12lf %d\n", indices[i], value);

        value *= 1.0001;
        value = DOUBLE_max(value, 0.0001);

        cut->indices[i] = indices[i];
        cut->pi[i] = - value;
    }

    cut->pi_zero = -1.0;

CLEANUP:
    LP_free(&lp);
    return rval;
}


#ifndef TEST_SOURCE
int GREEDY_generate_cut(int nrows,
                        struct Row **rows,
                        const char *column_types,
                        struct Row *cut)
{
    int rval = 0;
    double *f = 0;
    long max_nrays = 0;

    f = (double *) malloc(nrows * sizeof(double));
    abort_if(!f, "could not allocate f");

    for (int i = 0; i < nrows; i++)
    {
        f[i] = frac(rows[i]->pi_zero);
        if (DOUBLE_eq(f[i], 1.0)) f[i] = 0.0;

        max_nrays += rows[i]->nz;
    }

    int nz;
    int nrays;
    int *variable_to_ray = 0;
    int *indices = 0;
    double *extracted_rays = 0;
    double *scale = 0;

    double *beta = 0;
    double *scaled_rays = 0;

    int lfree_nrays = 0;
    double *lfree_rays = 0;

    variable_to_ray = (int *) malloc(max_nrays * sizeof(int));
    indices = (int *) malloc(max_nrays * sizeof(int));
    scale = (double *) malloc(max_nrays * sizeof(double));
    extracted_rays = (double *) malloc(nrows * max_nrays * sizeof(double));
    abort_if(!variable_to_ray, "could not allocate variable_to_ray");
    abort_if(!indices, "could not allocate indices");
    abort_if(!scale, "could not allocate scale");
    abort_if(!extracted_rays, "could not allocate extracted_rays");

    lfree_rays = (double*) malloc(nrows * (max_nrays + 100) * sizeof(double));
    abort_if(!lfree_rays, "could not allocate lfree_rays");

    rval = CG_extract_rays_from_rows(nrows, rows, extracted_rays, &nrays,
            variable_to_ray, scale, indices, &nz);
    abort_if(rval, "CG_extract_rays_from_rows failed");

    for (double norm_cutoff = 0.00; norm_cutoff <= 5.0; norm_cutoff+= 0.1)
    {
        lfree_nrays = 0;

        for (int i = 0; i < nrays; i++)
        {
            int keep = 1;
            double *r = &extracted_rays[nrows * i];

            for (int j = 0; j < lfree_nrays; j++)
            {
                double *q = &extracted_rays[nrows * j];
                double norm = 0;

                for(int k = 0; k < nrows; k++)
                    norm += fabs(r[k] - q[k]);
                
                if(norm <= norm_cutoff)
                {
                    keep = 0;
                    break;
                }
            }

            if(keep)
            {
                memcpy(&lfree_rays[nrows * lfree_nrays], r, nrows * sizeof(double));
                lfree_nrays++;
            }
        }

        log_debug("  norm_cutoff=%8.2lf nrays=%8d\n", norm_cutoff, lfree_nrays);

        if(lfree_nrays < MAX_N_RAYS) break;
    }

    if(ENABLE_LIFTING)
    {
        abort_if(nrows > 3, "not implemented");

        int n_extra_rays;
        double extra_rays[100];
        
        if(nrows == 2)
        {
            extra_rays[0] = 0.0001;
            extra_rays[1] = 0.0000;

            extra_rays[2] = -0.0001;
            extra_rays[3] =  0.0001;

            extra_rays[4] = -0.0001;
            extra_rays[5] = -0.0001;

            n_extra_rays = 3;
        }
        else if(nrows == 3)
        {
            extra_rays[0] = 0.0000;
            extra_rays[1] = 0.0000;
            extra_rays[2] = 0.0001;

            extra_rays[3] =  0.0001;
            extra_rays[4] =  0.0000;
            extra_rays[5] = -0.0001;
            
            extra_rays[6] = -0.0001;
            extra_rays[7] =  0.0000;
            extra_rays[8] = -0.0001;
            
            extra_rays[ 9] = -0.0001;
            extra_rays[10] = -0.0001;
            extra_rays[11] = -0.0001;

            n_extra_rays = 4;
        }


        for(int i = 0; i < n_extra_rays; i++)
        {
            double *r = &extra_rays[nrows * i];

            double scale;
            int found, index;

            rval = CG_find_ray(nrows, lfree_rays, lfree_nrays, r, &found,
                    &scale, &index);
            abort_if(rval, "CG_find_ray failed");

            if(!found) {
                memcpy(&lfree_rays[nrows * lfree_nrays], r, nrows *
                        sizeof(double));
                lfree_nrays++;
            }
        }
    }


    if(lfree_nrays < 3)
    {
        rval = ERR_NO_CUT;
        cut->pi = 0;
        cut->indices = 0;
        goto CLEANUP;
    }

    log_debug("Extracted %d rays\n", lfree_nrays);
    if_verbose_level
    {
        print_rays(variable_to_ray, indices, f, extracted_rays, scale, nrays,
                nz, nrows);
    }

    beta = (double *) malloc((lfree_nrays) * sizeof(double));
    abort_if(!beta, "could not allocate beta");

    log_verbose("Computing lattice-free set...\n");

    if(nrows == 2)
    {
        rval = sort_rays_angle(lfree_rays, lfree_nrays, beta);
        abort_if(rval, "sort_rays_angle failed");

        rval = GREEDY_2D_generate_cut(lfree_rays, lfree_nrays, f, beta);
        if(rval)
        {
            rval = ERR_NO_CUT;
            goto CLEANUP;
        }
    }
    else
    {
        rval = GREEDY_ND_generate_cut(nrows, lfree_nrays, f, lfree_rays, beta);
        if(rval)
        {
            rval = ERR_NO_CUT;
            goto CLEANUP;
        }
    }


    if(DUMP_CUT)
    {
        char filename[100];
        sprintf(filename, "cut-%03d.sage", DUMP_CUT_N++);

        time_printf("Writing %s...\n", filename);
        rval = GREEDY_write_sage_file(nrows, lfree_nrays, f, lfree_rays, beta,
                filename);
        abort_if(rval, "GREEDY_write_sage_file failed");
    }

    rval = create_cut(nrows, column_types, nrays, f, lfree_nrays, lfree_rays,
            beta, extracted_rays, nz, variable_to_ray, indices, scale, cut);
    abort_if(rval, "create_cut failed");

    if_verbose_level print_cut(cut);

CLEANUP:
    if (f) free(f);
    if (variable_to_ray) free(variable_to_ray);
    if (beta) free(beta);
    if (indices) free(indices);
    if (scale) free(scale);
    if (scaled_rays) free(scaled_rays);
    if (lfree_rays) free(lfree_rays);
    return rval;
}
#endif // TEST_SOURCE
