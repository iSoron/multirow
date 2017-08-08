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
#include <infinity/infinity-nd.h>
#include <infinity/infinity-2d.h>

/**
 * Auxiliary structure for sorting with qsort.
 */
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
 * @param rays the list to be sorted
 * @return zero if successful, non-zero otherwise
 */
static int sort_rays_by_angle(struct RayList *rays)
{
    int rval = 0;
    int nrays = rays->nrays;
    double *rays_copy = 0;
    struct SortPair *pairs = 0;

    pairs = (struct SortPair *) malloc(nrays * sizeof(struct SortPair));
    rays_copy = (double *) malloc(2 * nrays * sizeof(double));
    abort_if(!pairs, "could not allocate pairs");
    abort_if(!rays_copy, "could not allocate rays_copy");

    memcpy(rays_copy, rays->values, 2 * nrays * sizeof(double));

    for(int i = 0; i < nrays; i++)
    {
        pairs[i].index = i;
        pairs[i].data = &rays->values[2 * i];
    }

    qsort(pairs, (size_t) nrays, sizeof(struct SortPair),
            _qsort_cmp_rays_angle);

    rays->nrays = 0;
    for(int i = 0; i < nrays; i++)
        LFREE_push_ray(rays, &rays_copy[2 * pairs[i].index]);

CLEANUP:
    free(pairs);
    free(rays_copy);
    return rval;
}

/**
 * Creates an intersection cut from the given lattice-free set.
 *
 * @param tableau the tableau that was used to generate the model
 * @param map the mapping between the tableau and model
 * @param model the multi-row model that was used to generate the cut
 * @param lfree the lattice-free set
 * @param[out] cut the generated cut
 * @return zero if successful, non-zero otherwise
 */
static int create_cut_from_lfree(const struct Tableau *tableau,
                                const struct TableauModelMap *map,
                                const struct MultiRowModel *model,
                                const struct ConvLFreeSet *lfree,
                                struct Row *cut)
{
    int rval = 0;
    double *ray = 0;

    struct LP lp;
    int nvars = map->nvars;
    int nrows = tableau->nrows;
    const struct RayList *rays = &model->rays;

    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    memcpy(lfree->f, model->f, nrows * sizeof(double));

    rval = INFINITY_create_psi_lp(lfree, &lp);
    abort_if(rval, "create_psi_lp failed");

    ray = (double*) malloc(nrows * sizeof(double));
    abort_if(!ray, "could not allocate ray");

    cut->nz = nvars;
    for(int i = 0; i < nvars; i++)
    {
        double value, norm = 0;
        double *original_ray = LFREE_get_ray(rays, map->variable_to_ray[i]);
        char type = tableau->column_types[map->indices[i]];

        for (int j = 0; j < nrows; j++)
        {
            ray[j] = original_ray[j];
            norm += fabs(ray[j]);
        }

        if (norm < 0.001)
            for (int j = 0; j < nrows; j++)
                ray[j] *= 0.001 / norm;

        if(ENABLE_LIFTING && type == MILP_INTEGER)
        {
            rval = INFINITY_pi(nrows, ray, map->ray_scale[i], &lp, &value);
            abort_if(rval, "INFINITY_pi failed");
        }
        else
        {
            rval = INFINITY_psi(nrows, ray, map->ray_scale[i], &lp, &value);
            abort_if(rval, "INFINITY_psi failed");
        }

        log_verbose("   psi[%4d] = %20.12lf %d\n", map->indices[i], value);

        if_debug_level
        {
            time_printf("    q=[");
            for(int j = 0; j < nrows; j++)
                printf("%25.20lf ", ray[j] * map->ray_scale[i]);
            printf("] value=%25.20lf\n", value);
        }

        value = fmax(value, 1 / 1024.0);

//        value *= 1.001;
//        value = DOUBLE_max(value, 0.001);

        cut->indices[i] = map->indices[i];
        cut->pi[i] = -value;
    }

    cut->pi_zero = -1.0;

CLEANUP:
    if(ray) free(ray);
    LP_free(&lp);
    return rval;
}

/**
 * Given an original multi-row model, returns a simplified model with fewer
 * rays. For each subset of rays that are very close to each other, only one ray
 * is selected for the new model.
 *
 * @param original_model the original model
 * @param[out] filtered_model the simplified model
 * @return
 */
static int filter_model(const struct MultiRowModel *original_model,
                        struct MultiRowModel *filtered_model)
{
    int rval = 0;
    double *r = 0;
    double *f = filtered_model->f;
    int nrows = original_model->nrows;
    struct RayList *filtered_rays = &filtered_model->rays;
    const struct RayList *original_rays = &original_model->rays;

    memcpy(f, original_model->f, 2 * sizeof(double));
    for(int i = 0; i < nrows; i++)
    {
        f[i] = (ceil(f[i] * 128) / 128);
        if(f[i] <= 0.01) f[i] = 0;
    }

    r = (double*) malloc(nrows * sizeof(double));
    abort_if(!r, "could not allocate r");

    for(double norm_cutoff = 0.00; norm_cutoff <= 5.0; norm_cutoff += 0.1)
    {
        filtered_rays->nrays = 0;

        for(int i = 0; i < original_rays->nrays; i++)
        {
            int keep = 1;
            double norm = 0;

            memcpy(r, LFREE_get_ray(original_rays, i), nrows * sizeof(double));

            for(int j = 0; j < nrows; j++)
            {
                r[j] = (ceil(r[j] * 128) / 128);
                norm += fabs(r[j]);
            }

            if(DOUBLE_iszero(norm)) continue;

            for(int j = 0; j < (filtered_rays->nrays); j++)
            {
                double *q = LFREE_get_ray(filtered_rays, j);

                norm = 0;
                for(int k = 0; k < nrows; k++)
                    norm += fabs(r[k] - q[k]);

                if(norm <= norm_cutoff)
                {
                    keep = 0;
                    break;
                }
            }

            if(keep) LFREE_push_ray(filtered_rays, r);
        }

        log_verbose("  norm_cutoff=%8.2lf nrays=%8d\n", norm_cutoff,
                filtered_model->rays.nrays);

        if(filtered_rays->nrays < MAX_N_RAYS) break;
    }

CLEANUP:
    if(r) free(r);
    return rval;
}

/**
 * Appends a few extra rays to the model, so that the resulting
 * lattice-free set is a little larger than it would be otherwise.
 * Although this does not change the coefficients for the continuous
 * variables, it might help with lifting.
 *
 * The model is modified in-place.
 *
 * @param model the model to be modified
 * @return zero if successful, non-zero otherwise
 */
static int append_extra_rays(struct MultiRowModel *model)
{
    int rval = 0;
    int nrows = model->nrows;

    double *r = 0;
    double e = 1 / 128.0;

    r = (double*) malloc(nrows * sizeof(double));
    abort_if(!r, "could not allocate r");

    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < nrows; j++) r[j] = (i == j ? e : 0);
        LFREE_push_ray(&model->rays, r);

        for(int j = 0; j < nrows; j++) r[j] = (i == j ? -e : 0);
        LFREE_push_ray(&model->rays, r);
    }

    if(nrows == 2)
    {
        rval = sort_rays_by_angle(&model->rays);
        abort_if(rval, "sort_rays_by_angle failed");
    }

CLEANUP:
    if(r) free(r);
    return rval;
}

/**
 * Writes a given lattice-free set to a sage file.
 *
 * This file can be rendered by the script 'benchmark/scripts/render.sage'
 *
 * @param lfree the lattice-free set to be written
 * @param filename the name of the file
 * @return zero if successful, non-zero otherwise
 */
static int write_lfree_to_sage_file(const struct ConvLFreeSet *lfree,
                                    const char *filename)
{
    int rval = 0;
    int nrows = lfree->nrows;

    FILE *fsage = fopen(filename, "w");
    abort_iff(!fsage, "could not open %s", filename);

    fprintf(fsage, "f=vector([");
    for(int i = 0; i < nrows; i++)
        fprintf(fsage, "%.20lf,", lfree->f[i]);
    fprintf(fsage, "])\n");

    fprintf(fsage, "R=matrix([\n");
    for(int i = 0; i < lfree->rays.nrays; i++)
    {
        double *r = LFREE_get_ray(&lfree->rays, i);
        fprintf(fsage, "    [");
        for(int j = 0; j < nrows; j++)
            fprintf(fsage, "%.20lf,", r[j]);
        fprintf(fsage, "],\n");
    }
    fprintf(fsage, "])\n");

    fprintf(fsage, "pi=vector([\n");
    for(int k = 0; k < lfree->rays.nrays; k++)
        fprintf(fsage, "    %.12lf,\n", 1 / lfree->beta[k]);
    fprintf(fsage, "])\n");

CLEANUP:
    fclose(fsage);
    return rval;
}

static int dump_cut(const struct ConvLFreeSet *lfree)
{
    int rval = 0;

    char filename[100];
    sprintf(filename, "cut-%03d.sage", DUMP_CUT_N++);

    time_printf("Writing %s...\n", filename);
    rval = write_lfree_to_sage_file(lfree, filename);
    abort_if(rval, "write_lfree_to_sage_file failed");

CLEANUP:
    return rval;
}

/**
 * Given a tableau, extracts the original multi-row model, along with
 * its mapping, in addition to a simplified version of this model.
 *
 * @param tableau the tableau
 * @param original_model the original multi-row model obtained from the tableau
 * @param filtered_model the simplified multi-row model obtained from the tableau
 * @param original_map the mapping between the tableau and the original model
 * @return zero if successful, non-zero otherwise
 */
static int extract_models(const struct Tableau *tableau,
                          struct MultiRowModel *original_model,
                          struct MultiRowModel *filtered_model,
                          struct TableauModelMap *original_map)
{
    int rval = 0;

    rval = CG_extract_f_from_tableau(tableau, original_model->f);
    abort_if(rval, "CG_extract_f_from_tableau failed");

    rval = CG_extract_model(tableau, original_map, original_model);
    abort_if(rval, "CG_extract_rays_from_rows failed");

    rval = filter_model(original_model, filtered_model);
    abort_if(rval, "filter_model failed");

    if(tableau->nrows == 2)
    {
        rval = sort_rays_by_angle(&filtered_model->rays);
        abort_if(rval, "sort_rays_by_angle failed");
    }

CLEANUP:
    return rval;
}

#ifndef TEST_SOURCE

/**
 * Generates the infinity cut for a given tableau.
 *
 * @param tableau the tableau that should be used to generate the cut
 * @param[out] cut the resulting infinity cut
 * @return zero if successful, non-zero otherwise
 */
int INFINITY_generate_cut(const struct Tableau *tableau, struct Row *cut)
{
    int rval = 0;
    int max_nrays = CG_total_nz(tableau) + 100;

    struct MultiRowModel original_model;
    struct MultiRowModel filtered_model;
    struct TableauModelMap original_map;
    struct ConvLFreeSet lfree;

    rval = CG_init_model(&original_model, tableau->nrows, max_nrays);
    abort_if(rval, "CG_init_model failed");

    rval = CG_init_model(&filtered_model, tableau->nrows, max_nrays);
    abort_if(rval, "CG_init_model failed");

    rval = CG_init_map(&original_map, max_nrays, tableau->nrows);
    abort_if(rval, "CG_init_map failed");

    rval = extract_models(tableau, &original_model, &filtered_model,
            &original_map);
    abort_if(rval, "extract_models failed");

    if(filtered_model.rays.nrays < 3)
    {
        rval = ERR_NO_CUT;
        goto CLEANUP;
    }

    rval = LFREE_init_conv(&lfree, tableau->nrows, max_nrays);
    abort_if(rval, "LFREE_init_conv failed");

//    if(tableau->nrows == 2)
//        rval = INFINITY_2D_generate_lfree(&filtered_model, &lfree);
//    else
//        rval = INFINITY_ND_generate_lfree(&filtered_model, &lfree);
//
//    if(rval)
//    {
//        rval = ERR_NO_CUT;
//        goto CLEANUP;
//    }
//
//    for(int i = 0; i < filtered_model.rays.nrays; i++)
//    {
//        double *r = LFREE_get_ray(&filtered_model.rays, i);
//        for(int j = 0; j < tableau->nrows; j++)
//            r[j] *= lfree.beta[j];
//    }

    rval = append_extra_rays(&filtered_model);
    abort_if(rval, "append_extra_rays failed");

//    if(tableau->nrows == 2)
//        rval = INFINITY_2D_generate_lfree(&filtered_model, &lfree);
//    else
        rval = INFINITY_ND_generate_lfree(&filtered_model, &lfree);

    if(rval)
    {
        rval = ERR_NO_CUT;
        goto CLEANUP;
    }

    if(SHOULD_DUMP_CUTS)
    {
        rval = dump_cut(&lfree);
        abort_if(rval, "dump_cut failed");
    }

    rval = create_cut_from_lfree(tableau, &original_map, &original_model,
            &lfree, cut);
    abort_if(rval, "create_cut_from_lfree failed");

CLEANUP:
    LFREE_free_conv(&lfree);
    CG_free_model(&original_model);
    CG_free_model(&filtered_model);
    CG_free_map(&original_map);
    return rval;
}

#endif // TEST_SOURCE