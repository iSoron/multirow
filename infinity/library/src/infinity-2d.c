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
#include <stdio.h>
#include <stdlib.h>

#include <multirow/geometry.h>
#include <multirow/double.h>
#include <multirow/util.h>
#include <multirow/cg.h>

#include <infinity/infinity-2d.h>

static int get_bounding_box(int nrows,
                            int nrays,
                            const double *rays,
                            const double *bounds,
                            double epsilon,
                            double *lb,
                            double *ub)
{
    for(int i = 0; i < nrows; i++)
    {
        ub[i] = 0;
        lb[i] = 0;
    }

    for (int i=0; i < nrays; i++)
    {
        double e = fmin(epsilon, bounds[i]);

        for(int j = 0; j < nrows; j++)
        {
            double rij = rays[nrows * i + j] * e;
            lb[j] = fmin(lb[j], floor(rij));
            ub[j] = fmax(ub[j], ceil(rij));
        }
    }

    return 0;
}

struct LatticeSequence
{
    int square;
    int direction;
    int steps;
    int eol;

    int i;
    int j;
};

static void lattice_sequence_init(struct LatticeSequence *seq)
{
    seq->steps = 0;
    seq->square = 0;
    seq->direction = 0;
    seq->i = 0;
    seq->j = 0;
    seq->eol = 0;
}

static int next_lattice_point(struct LatticeSequence *seq,
                              const double *lb,
                              const double *ub)
{
    int rval = 0;

    seq->eol = 1;

    for(int k = 0; k < 8 * (seq->square+1); k++)
    {
        if(seq->steps > 0)
        {
            seq->steps--;

            switch(seq->direction)
            {
                case 3:
                    seq->i++;
                    break;

                case 2:
                    seq->j--;
                    break;

                case 1:
                    seq->i--;
                    break;

                case 0:
                    seq->j++;
                    break;
            }
        }

        if(seq->steps == 0)
        {
            if(seq->direction > 0)
            {
                seq->direction--;
                seq->steps = 2 * seq->square;
            }
            else
            {
                seq->square++;
                seq->direction = 3;
                seq->steps = 2 * seq->square;
                seq->i--;
                seq->j++;
            }
        }

        if(seq->i >= lb[0] && seq->i <= ub[0] && seq->j >= lb[1] &&
                seq->j <= ub[1])
        {
            seq->eol = 0;
            break;
        }
    }

    CLEANUP:
    return rval;
}



/*
 * Returns lambda1 and lambda2 such that p = r1 * lambda1 + r2 * lambda2
 */
static void find_linear_combination(const double *r1,
                                    const double *r2,
                                    const double *p,
                                    double *lambda1,
                                    double *lambda2)
{
    int rval = 0;

    double den = (r1[0] * r2[1] - r2[0] * r1[1]);

    if(DOUBLE_iszero(den)) den = 0.0;

    *lambda1 = (r2[1] * p[0] - r2[0] * p[1]) / den;
    *lambda2 = (-r1[1] * p[0] + r1[0] * p[1]) / den;
}


static int generate_split(const double *f,
                          const double *d,
                          double *pi,
                          double *pi_zero,
                          long max_den)
{
    int rval = 0;

    Rational d1, d2;
    double m;

    rval = DOUBLE_to_rational(d[0], 10, d1);
    abort_if(rval, "DOUBLE_to_rational failed");

    rval = DOUBLE_to_rational(d[1], 10, d2);
    abort_if(rval, "DOUBLE_to_rational failed");

    m = lcm(d1->den, d2->den);

    d1->den = m / d1->den;
    d1->num *= d1->den;

    d2->den = m / d2->den;
    d2->num *= d2->den;

    m = gcd(d1->num, d2->num);
    if(m != 0)
    {
        d1->num /= m;
        d2->num /= m;

        pi[0] = (double) d2->num;
        pi[1] = - ((double) d1->num);
        *pi_zero = floor(pi[0] * f[0] + pi[1] * f[1]);
    }
    else
    {
        pi[0] = pi[1] = INFINITY;
        *pi_zero = INFINITY;
    }

    CLEANUP:
    return rval;
}

/*
 * Receives a list of rays r1,...,rn and a point p. Returns i1,i2 such that
 * p belongs to cone(r_i1,r_i2). Also returns lambda1, lambda2 such that
 * p = r_i1 * lambda1 + r_i2 * lambda2.
 *
 * The rays must be sorted in clockwise order.
 */
static int find_containing_cone(const double *rays,
                                const int nrays,
                                const double *p,
                                int *index1,
                                int *index2,
                                double *lambda1,
                                double *lambda2)
{
    int rval = 0;

    int i1, i2;

    for (i1 = 0; i1 < nrays; i1++)
    {
        i2 = (i1 + 1) % nrays;

        const double *r1 = &rays[i1 * 2];
        const double *r2 = &rays[i2 * 2];

        double at1 = atan2(r1[1], r1[0]);
        double at2 = atan2(r2[1], r2[0]) - at1;
        if (at2 > 0) at2 -= 2 * M_PI;
        if (at2 <= - M_PI)
        {
            log_verbose("    discarding obtuse rays\n");
            log_verbose("      r1=%.12lf %.12lf\n", r1[0], r1[1]);
            log_verbose("      r2=%.12lf %.12lf\n", r2[0], r2[1]);
            continue;
        }

        find_linear_combination(r1, r2, p, lambda1, lambda2);

        log_verbose("      r1=%.12lf %.12lf\n", r1[0], r1[1]);
        log_verbose("      r2=%.12lf %.12lf\n", r2[0], r2[1]);
        log_verbose("      %.8lf %.8lf\n", *lambda1, *lambda2);

        if(DOUBLE_iszero(*lambda1)) *lambda1 = 0.0;
        if(DOUBLE_iszero(*lambda2)) *lambda2 = 0.0;

        if ((*lambda1) >= 0 && (*lambda2) >= 0 && (*lambda1) <= 1e9 &&
            (*lambda2) <= 1e9)
        {
            log_verbose("        match!\n");
            break;
        }
    }

    if (i1 == nrays)
        i1 = i2 = -1;

    *index1 = i1;
    *index2 = i2;

CLEANUP:
    return rval;
}


/*
 * Find lambda such that p lies on the line connecting lambda * r1 and lambda * r2
 */
static int scale_cone_to_point(const double *r1,
                               const double *r2,
                               const double *p,
                               double *lambda)
{
    int rval = 0;

    double a = r1[0], b = r1[1];
    double c = r2[0], d = r2[1];
    double den = (b * c - a * d);

    //abort_iff(fabs(den) < 1e-9, "rays cannot be parallel (den=%.12lf)", den);

    *lambda = p[0] * (b - d) - p[1] * (a - c);
    *lambda /= den;

CLEANUP:
    return rval;
}

/*
 * Find lambda such that p lies in the line connecting r1 and lambda * r2
 */
static int shear_cone_to_point(const double *r1,
                               const double *r2,
                               const double *p,
                               double *lambda)
{
    int rval = 0;

    double a = r1[0], b = r1[1];
    double c = r2[0], d = r2[1];
    double den = d * (p[0] - a) - c * (p[1] - b);

    *lambda = b * p[0] - a * p[1];
    *lambda /= den;

CLEANUP:
    return rval;
}

static int scale_to_chull(double *rays, int nrays, double *scale)
{
    int rval = 0;

    double *rays_extended = 0;
    double *vertices = 0;
    int nvertices;

    rays_extended = (double*) malloc(2 * (nrays + 1) * sizeof(double));
    vertices = (double*) malloc(2 * (nrays + 1) * sizeof(double));
    abort_if(!rays_extended, "could not allocate rays_extended");
    abort_if(!vertices, "could not allocate vertices");

    memcpy(rays_extended, rays, 2 * nrays * sizeof(double));
    rays_extended[2*nrays] = 0;
    rays_extended[2*nrays + 1] = 0;

    rval = chull_2d(rays_extended, nrays + 1, vertices, &nvertices);
    abort_if(rval, "chull_2d failed");

    for(int i = 0; i < nrays; i++)
        scale[i] = 1.0;

    log_verbose("  convex hull:\n");
    for(int i = 0; i < nvertices; i++)
    {
        log_verbose("    v%-3d: %20.8lf %20.8lf\n", i, vertices[2 * i],
                vertices[2 * i + 1]);
    }

    log_verbose("  rays:\n");
    for(int i = 0; nrays >= 3 && i < nrays; i++)
    {
        int i1, i2;
        double lambda1, lambda2, mu;

        rval = find_containing_cone(vertices, nvertices, &rays[2*i], &i1, &i2,
                &lambda1, &lambda2);
        abort_if(rval, "find_containing_cone failed");
        log_verbose("%.8lf %.8lf\n", lambda1, lambda2);

        if(i1 < 0 || i2 < 0) continue;

        if(DOUBLE_iszero(lambda1))
        {
            mu = lambda2;
        }
        else if(DOUBLE_iszero(lambda2))
        {
            mu = lambda1;
        }
        else
        {
            rval = scale_vector_to_line(&vertices[2*i1], &vertices[2*i2],
                    &rays[2*i], &mu);
            abort_if(rval, "scale_vector_to_line failed");
        }

        abort_if(!isfinite(mu), "mu should be finite");

        //log_verbose("  r%-3d: %.2lf %.2lf %.2lf\n", i, rays[2*i], rays[2*i+1], mu);

        rays[2*i] *= mu;
        rays[2*i+1] *= mu;
        scale[i] = mu;

        log_verbose("    r%-3d: %20.12lf %20.12lf (scale %.8lf)\n", i, rays[2*i],
                rays[2*i+1], scale[i]);
    }

CLEANUP:
    if(rays_extended) free(rays_extended);
    if(vertices) free(vertices);
    return rval;
}

static int bound(const double *rays,
                 const double *bounds,
                 int nrays,
                 const double *f,
                 const double *p,
                 double *epsilon,
                 double *v1,
                 double *v2,
                 int *index1,
                 int *index2)
{
    int rval = 0;

    int i1, i2, iexact = -1;

    double e1, e2;
    const double *r1, *r2;
    double lambda1, lambda2;

    double pp[2] = { p[0] - f[0], p[1] - f[1] };

    rval = find_containing_cone(rays, nrays, pp, &i1, &i2, &lambda1, &lambda2);
    abort_if(rval, "find_containing_cone failed");

    if(i1 < 0 || i2 < 0)
    {
        log_verbose("    no cone\n");
        *epsilon = INFINITY;
        goto CLEANUP;
    }

    if(DOUBLE_iszero(lambda1)) iexact = i2;
    if(DOUBLE_iszero(lambda2)) iexact = i1;

    if(iexact >= 0)
    {
        int inext = (iexact + 1) % nrays;
        int iprev = (iexact + (nrays-1)) % nrays;

        double mu1, mu2;
        find_linear_combination(&rays[2 * iprev], &rays[2 * inext],
                &rays[2 * iexact], &mu1, &mu2);

        log_verbose("    mu=%.12lf %.12lf\n", mu1, mu2);

        int should_enlarge_cone = 1;
        if(!DOUBLE_leq(bounds[iexact], lambda1+lambda2)) should_enlarge_cone = 0;
        if(!DOUBLE_geq(mu1, 0) || !DOUBLE_geq(mu2, 0)) should_enlarge_cone = 0;
        if(isinf(mu1) || isinf(mu2)) should_enlarge_cone = 0;
        if(fabs(mu1) > 1e9 || fabs(mu2) > 1e9) should_enlarge_cone = 0;

        if(should_enlarge_cone)
        {
            i1 = iprev;
            i2 = inext;
        }
        else
        {
            v1[0] = v1[1] = INFINITY;
            v2[0] = v2[1] = INFINITY;
            *epsilon = (lambda1 + lambda2);
            if(DOUBLE_leq(bounds[iexact], *epsilon))
                *epsilon = INFINITY;

            *index1 = *index2 = iexact;
            goto CLEANUP;
        }

    }

    if(bounds[i1] > bounds[i2])
    {
        swap(i1, i2, int);
        swap(lambda1, lambda2, double);
    }

    *index1 = i1;
    *index2 = i2;

    r1 = &rays[2 * i1];
    r2 = &rays[2 * i2];

    log_verbose("      r%-4d  %20.12lf %20.12lf\n", i1, r1[0], r1[1]);
    log_verbose("      r%-4d  %20.12lf %20.12lf\n", i2, r2[0], r2[1]);
    log_verbose("      pp     %20.12lf %20.12lf\n", pp[0], pp[1]);
    log_verbose("      lambda %20.12lf %20.12lf\n", lambda1, lambda2);

    double r1bound[2] = { r1[0] * bounds[i1], r1[1] * bounds[i1] };

    rval = scale_cone_to_point(r1, r2, pp, &e1);
    abort_if(rval, "scale_cone_to_point failed");

    rval = shear_cone_to_point(r1bound, r2, pp, &e2);
    abort_if(rval, "scale_cone_to_point failed");

    log_verbose("    e1=%20.12lf\n", e1);
    log_verbose("    e2=%20.12lf\n", e2);
    log_verbose("    b1=%20.12lf\n", bounds[i1]);
    log_verbose("    b2=%20.12lf\n", bounds[i2]);

    switch(DOUBLE_cmp(e1, bounds[i1]))
    {
        case -1:
            *epsilon = e1;
            v1[0] = r1[0] * e1;
            v1[1] = r1[1] * e1;
            v2[0] = r2[0] * e1;
            v2[1] = r2[1] * e1;
            break;

        case 0:
        case 1:
            if(DOUBLE_geq(e2, bounds[i2]))
            {
                *epsilon = INFINITY;
            }
            else
            {
                *epsilon = e2;
                v1[0] = r1[0] * bounds[i1];
                v1[1] = r1[1] * bounds[i1];
                v2[0] = r2[0] * e2;
                v2[1] = r2[1] * e2;
            }
    }

CLEANUP:
    return rval;
}


#ifndef TEST_SOURCE


int INFINITY_2D_generate_cut(const struct MultiRowModel *model, double *bounds)
{
    log_verbose("INFINITY_2D_generate_cut\n");
    int rval = 0;
    int count = 0;
    int nrays = model->nrays;
    double *f = model->f;

    double *scale = 0;
    double *rays = 0;

    double lb[2], ub[2];

    for (int i = 0; i < nrays; i++)
        bounds[i] = GREEDY_BIG_E;

    scale = (double*) malloc(nrays * sizeof(double));
    rays = (double*) malloc(2 * nrays * sizeof(double));
    abort_if(!rays, "could not allocate rays");
    abort_if(!scale, "could not allocate scale");

    memcpy(rays, model->rays, 2 * nrays * sizeof(double));

    rval = scale_to_chull(rays, nrays, scale);
    abort_if(rval, "scale_to_chull failed");

    long seq_count = 0;

    while(1)
    {
        log_verbose("  starting iteration %d...\n", count);

        abort_if(count++ > 2 * nrays, "infinite loop");

        rval = get_bounding_box(2, nrays, rays, bounds, GREEDY_BIG_E, lb, ub);
        abort_if(rval, "get_bounding_box failed");

        log_verbose("    box=[%.2lf %.2lf] [%.2lf %.2lf]\n", lb[0], ub[0], lb[1], ub[1]);

        int best_i1, best_i2;
        double best_v1[2];
        double best_v2[2];
        double best_p[2];
        double best_epsilon = INFINITY;

        struct LatticeSequence seq;
        lattice_sequence_init(&seq);

        while(!seq.eol)
        {
            seq_count++;
            double p[2] = { seq.i, seq.j };
            double v1[2], v2[2], epsilon;
            int i1, i2;

            log_verbose("    p=%.2lf %.2lf\n", p[0], p[1]);

            rval = bound(rays, bounds, nrays, f, p, &epsilon, v1, v2, &i1, &i2);
            abort_if(rval, "bound failed");

            log_verbose("     epsilon=%.2lf\n", epsilon);

            if(epsilon >= 0 && epsilon < best_epsilon)
            {
                log_verbose("    found smaller epsilon: %.8lf\n", epsilon);

                rval = get_bounding_box(2, nrays, rays, bounds, epsilon, lb, ub);
                abort_if(rval, "get_bounding_box failed");

                log_verbose("      p=%.2lf %.2lf\n", p[0], p[1]);
                log_verbose("      box=[%.2lf %.2lf] [%.2lf %.2lf]\n", lb[0], ub[0],
                        lb[1], ub[1]);
                log_verbose("      v1=%12.8lf %12.8lf\n", v1[0], v1[1]);
                log_verbose("      v2=%12.8lf %12.8lf\n", v2[0], v2[1]);
                log_verbose("      rays=[%d %d]\n", i1, i2);

                best_epsilon = epsilon;
                best_v1[0] = v1[0];
                best_v1[1] = v1[1];
                best_v2[0] = v2[0];
                best_v2[1] = v2[1];
                best_p[0] = p[0];
                best_p[1] = p[1];
                best_i1 = i1;
                best_i2 = i2;
            }

            next_lattice_point(&seq, lb, ub);
            if(seq_count > MAX_LATTICE_POINTS)
            {
                rval = ERR_MIP_TIMEOUT;
                goto CLEANUP;
            }
        }

        if(isinf(best_epsilon))
        {
            log_verbose("    best_epsilon is infinity\n");
            break;
        }

        log_verbose("  updating bounds\n");
        if(isinf(best_v1[0]))
        {
            bounds[best_i1] = best_epsilon;
            log_verbose("    bounds[%d]=%.8lf (exact)\n", best_i1, best_epsilon);
        }
        else
        {
            log_verbose("    v1=%.8lf %.8lf\n", best_v1[0], best_v1[1]);
            log_verbose("    v2=%.8lf %.8lf\n", best_v2[0], best_v2[1]);
            log_verbose("    i=%d %d\n", best_i1, best_i2);

            for (int i = 0; i < nrays; i++)
            {
                double lambda;

                rval = scale_vector_to_line(best_v1, best_v2, &rays[2 * i], &lambda);
                abort_if(rval, "scale_vector_to_line failed");

                if(!DOUBLE_geq(lambda, 0)) continue;

                bounds[i] = fmin(bounds[i], lambda);
                log_verbose("    bounds[%d]=%.8lf\n", i, bounds[i]);
            }
        }

        //if(count > 0)
        //{
        //    for (int i = 0; i < nrays; i++)
        //       bounds[i] = fmin(bounds[i], best_epsilon);

        //    break;
        //}

        int is_split;

        for (int k = 0; k < nrays; k++)
        {
            if(bounds[k] < 100) continue;
            is_split = 1;

            double *split_direction = &rays[2 * k];

            log_verbose("  split_direction=%.2lf %.2lf\n", split_direction[0],
                    split_direction[1]);

            double pi[2];
            double pi_zero;

            rval = generate_split(f, split_direction, pi, &pi_zero, 10);
            abort_if(rval, "generate_split failed");

            log_verbose("    pi=%.2lf %.2lf\n", pi[0], pi[1]);
            log_verbose("    pi_zero=%.2lf\n", pi_zero);

            if(isinf(pi[0]))
            {
                is_split = 0;
                break;
            }

            double lhs;

            // reject splits that have f on the boundary
            lhs = f[0] * pi[0] + f[1] * pi[1];
            if(DOUBLE_eq(pi_zero, lhs) || DOUBLE_eq(lhs, pi_zero+1))
            {
                log_verbose("    split rejected\n");
                is_split = 0;
            }

            for (int i = 0; i < nrays && is_split; i++)
            {
                const double *r = &rays[2 * i];

                lhs = (f[0] + r[0] * bounds[i]) * pi[0];
                lhs += (f[1] + r[1] * bounds[i]) * pi[1];

                if (!(DOUBLE_leq(pi_zero, lhs) && DOUBLE_leq(lhs, pi_zero+1)))
                {
                    log_verbose("    point %.4lf %.4lf falls outside of the split\n",
                            f[0] + r[0]*bounds[i], f[1] + r[1] * bounds[i]);
                    is_split = 0;
                }
            }

            if(is_split)
            {
                log_verbose("  split confirmed. stopping.\n");
                break;
            }
        }

        if(is_split) break;
    }

    for(int i=0; i<nrays; i++)
        bounds[i] *= scale[i];

    CLEANUP:
    if(scale) free(scale);
    if(rays) free(rays);
    return rval;
}

#endif // TEST_SOURCE
