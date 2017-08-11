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

#include <gtest/gtest.h>

#define TEST_SOURCE

extern "C" {
#include <multirow/lp.h>
#include <multirow/util.h>
#include <multirow/lfree2d.h>
#include <multirow/cg.h>
#include <infinity/infinity-nd.h>
#include "../src/infinity-nd.c"
}


TEST(InfinityNDTest, find_violated_cone_test)
{
    int rval = 0;

    int nrows = 2;
    int nrays = 4;
    double f[] = { 0.5, 0.5 };
    double rays[] =
    {
        -1.0,  1.0,
         1.0,  1.0,
         1.0, -1.0,
        -1.0, -1.0
    };

    double x[] = { 1.0, 0.5 };
    double beta[] = { INFINITY, INFINITY, INFINITY, INFINITY };

    double sbar[nrays];
    int rx[nrays];
    int found;

    rval = find_violated_cone(nrows, nrays, f, rays, x, beta, 1.0, rx, sbar,
            &found);
    abort_if(rval, "find_violated_cone failed");

    EXPECT_TRUE(found);
    EXPECT_FALSE(rx[0]);
    EXPECT_TRUE(rx[1]);
    EXPECT_TRUE(rx[2]);
    EXPECT_FALSE(rx[3]);

    rval = find_violated_cone(nrows, nrays, f, rays, x, beta, 0.5, rx, sbar,
            &found);
    abort_if(rval, "find_violated_cone failed");

    EXPECT_FALSE(found);

    CLEANUP:
    if(rval) FAIL();
}

TEST(InfinityNDTest, find_tight_rays_test_1)
{
    int rval = 0;

    double f[] = { 0.5, 0.5 };
    double rays[] =
    {
         1.0,  1.0,
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
         0.0,  1.0,
         1.0,  0.0
    };
    double beta[] = { INFINITY, INFINITY, INFINITY, INFINITY, INFINITY,
                      INFINITY };
    double epsilon = 0.5;
    double x[] = { 1.0, 1.0 };

    int tx[6];

    rval = find_tight_rays(2, 6, f, rays, x, beta, epsilon, tx);
    abort_if(rval, "find_tight_rays failed");

    EXPECT_TRUE(tx[0]);
    EXPECT_FALSE(tx[1]);
    EXPECT_FALSE(tx[2]);
    EXPECT_FALSE(tx[3]);
    EXPECT_FALSE(tx[4]);
    EXPECT_FALSE(tx[5]);

    CLEANUP:
    if(rval) FAIL();
}

TEST(InfinityNDTest, find_tight_rays_test_2)
{
    int rval = 0;

    double f[] = { 0.5, 0.5 };
    double rays[] =
    {
         1.0,  1.0,
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
         0.0,  1.0,
         1.0,  0.0
    };
    double beta[] = { 0.5, 0.5, 0.5, 0.5, INFINITY, INFINITY };
    double epsilon = 1.0;
    double x[] = { 1.0, 1.0 };

    int tx[6];

    rval = find_tight_rays(2, 6, f, rays, x, beta, epsilon, tx);
    abort_if(rval, "find_tight_rays failed");

    EXPECT_TRUE(tx[0]);
    EXPECT_FALSE(tx[1]);
    EXPECT_FALSE(tx[2]);
    EXPECT_FALSE(tx[3]);
    EXPECT_TRUE(tx[4]);
    EXPECT_TRUE(tx[5]);

    CLEANUP:
    if(rval) FAIL();
}

TEST(InfinityNDTest, cone_bound_test_1)
{
    int rval = 0;

    double f[] = { 0.5, 0.5 };
    double rays[] =
    {
         1.0,  1.0,
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
         0.0,  1.0,
         1.0,  0.0
    };
    double beta[] = { INFINITY, INFINITY, INFINITY, INFINITY, INFINITY,
                      INFINITY };

    double x[] = { 1.0, 1.0 };
    int rx1[] = { 1, 0, 0, 0, 0, 0 };
    int rx2[] = { 0, 0, 0, 0, 1, 1 };

    double epsilon;

    rval = cone_bound(2, 6, f, rays, rx1, x, beta, &epsilon);
    abort_if(rval, "cone_bound failed");
    EXPECT_NEAR(0.5, epsilon, 1e-6);

    rval = cone_bound(2, 6, f, rays, rx2, x, beta, &epsilon);
    abort_if(rval, "cone_bound failed");
    EXPECT_NEAR(1.0, epsilon, 1e-6);

    CLEANUP:
    if(rval) FAIL();
}

TEST(InfinityNDTest, cone_bound_test_2)
{
    int rval = 0;

    double f[] = { 0.0, 0.0 };
    double rays[] =
    {
        0.0, 1.0,
        1.0, 0.0
    };
    int rx[] = { 1, 1 };

    double beta1[] = { 100, 100 };
    double beta2[] = { 0.5, 100 };
    double beta3[] = { 0.5, 1.0 };

    double x1[] = { 0.5, 0.5 };
    double x2[] = { 0.5, 0.25 };

    double epsilon;

    rval = cone_bound(2, 2, f, rays, rx, x1, beta1, &epsilon);
    abort_if(rval, "cone_bound failed");
    EXPECT_NEAR(1.0, epsilon, 1e-6);

    rval = cone_bound(2, 2, f, rays, rx, x1, beta2, &epsilon);
    abort_if(rval, "cone_bound failed");
    EXPECT_EQ(INFINITY, epsilon);

    rval = cone_bound(2, 2, f, rays, rx, x2, beta2, &epsilon);
    abort_if(rval, "cone_bound failed");
    EXPECT_NEAR(1.0, epsilon, 1e-6);

    rval = cone_bound(2, 2, f, rays, rx, x2, beta3, &epsilon);
    abort_if(rval, "cone_bound failed");
    EXPECT_EQ(INFINITY, epsilon);

CLEANUP:
    if(rval) FAIL();
}

TEST(InfinityNDTest, bound_test_1)
{
    int rval = 0;

    double f[] = { 0.5, 0.5 };
    double rays[] =
    {
         1.0,  1.0,
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
         0.0,  1.0,
         1.0,  0.0
    };
    double beta1[] = { INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY };
    double beta2[] = { 0.5, 0.5, 0.5, 0.5, INFINITY, INFINITY };
    double beta3[] = { 0.5, 0.5, 0.5, 0.5, 1.0, 1.0 };
    double x[] = { 1.0, 1.0 };

    double epsilon;
    int tx[6];

    rval = bound(2, 6, f, rays, x, beta1, &epsilon, tx);
    abort_if(rval, "bound failed");
    EXPECT_NEAR(epsilon, 0.5, 1e-6);
    EXPECT_TRUE(tx[0]);
    EXPECT_FALSE(tx[1]);
    EXPECT_FALSE(tx[2]);
    EXPECT_FALSE(tx[3]);
    EXPECT_FALSE(tx[4]);
    EXPECT_FALSE(tx[5]);

    rval = bound(2, 6, f, rays, x, beta2, &epsilon, tx);
    abort_if(rval, "bound failed");
    EXPECT_NEAR(epsilon, 1.0, 1e-6);
    EXPECT_TRUE(tx[0]);
    EXPECT_FALSE(tx[1]);
    EXPECT_FALSE(tx[2]);
    EXPECT_FALSE(tx[3]);
    EXPECT_TRUE(tx[4]);
    EXPECT_TRUE(tx[5]);

    rval = bound(2, 6, f, rays, x, beta3, &epsilon, tx);
    abort_if(rval, "bound failed");
    EXPECT_EQ(epsilon, INFINITY);
    EXPECT_FALSE(tx[0]);
    EXPECT_FALSE(tx[1]);
    EXPECT_FALSE(tx[2]);
    EXPECT_FALSE(tx[3]);
    EXPECT_FALSE(tx[4]);
    EXPECT_FALSE(tx[5]);

CLEANUP:
    if(rval) FAIL();
}

TEST(InfinityNDTest, psi_test)
{
    int rval = 0;

    double f[] = { 0.5, 0.5 };
    double rays[] =
    {
         1.0,  1.0,
         1.0, -1.0,
        -1.0, -1.0,
        -1.0,  1.0,
         0.0,  1.0,
         1.0,  0.0
    };
    double beta[] = { 0.5, 0.5, 0.5, 0.5, 1.0, 1.0 };

    double q1[] = { 1.0, 1.0 };
    double q2[] = { -2.0, 0.0 };

    struct ConvLFreeSet lfree;
    lfree.f = f;
    lfree.beta = beta;
    lfree.rays.nrays = 6;
    lfree.rays.values = rays;
    lfree.nrows = lfree.rays.dim = 2;

    double value;
    struct LP lp;

    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    rval = INFINITY_create_psi_lp(&lfree, &lp);
    abort_if(rval, "INFINITY_create_psi_lp failed");

    rval = INFINITY_psi(2, q1, 1.0, &lp, &value);
    abort_if(rval, "GREDDY_ND_psi failed");
    EXPECT_NEAR(value, 2.0, 1e-6);

    rval = INFINITY_psi(2, q2, 2.0, &lp, &value);
    abort_if(rval, "GREDDY_ND_psi failed");
    EXPECT_NEAR(value, 8.0, 1e-6);

CLEANUP:
    LP_free(&lp);
    if(rval) FAIL();
}

TEST(InfinityNDTest, psi_test_2)
{
    int rval = 0;

    double f[] = { 0.5, 0.5, 0.5 };
    double rays[] =
    {
        -0.5, -0.5, -0.5,
        -0.5, -0.5,  0.5,
        -0.5,  0.5, -0.5,
        -0.5,  0.5,  0.5,
         0.5, -0.5, -0.5,
         0.5, -0.5,  0.5,
         0.5,  0.5, -0.5,
         0.5,  0.5,  0.5,

    };

    double beta[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    double q1[] = { 0.5, 0.5, 0.5 };
    double q2[] = { 1, 0, 0 };

    struct ConvLFreeSet lfree;
    lfree.f = f;
    lfree.beta = beta;
    lfree.rays.nrays = 8;
    lfree.rays.values = rays;
    lfree.nrows = lfree.rays.dim = 3;

    double value;
    struct LP lp;

    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    rval = INFINITY_create_psi_lp(&lfree, &lp);
    abort_if(rval, "INFINITY_create_psi_lp failed");

    rval = INFINITY_psi(3, q1, 1.0, &lp, &value);
    abort_if(rval, "GREDDY_ND_psi failed");
    EXPECT_NEAR(value, 1.0, 1e-6);

    rval = INFINITY_psi(3, q2, 1.0, &lp, &value);
    abort_if(rval, "GREDDY_ND_psi failed");
    EXPECT_NEAR(value, 2.0, 1e-6);

CLEANUP:
    if(rval) FAIL();
}

TEST(InfinityNDTest, psi_test_3)
{
    int rval = 0;

    double f[] = { 0.671875, 0.671875 };
    double rays[] = {
            -0.007812500000,  0.000000000000,
            -0.039062500000,  0.046875000000,
             0.000000000000,  0.046875000000,
             0.046875000000,  0.000000000000,
             0.000000000000, -0.039062500000
    };

    double beta[] = {
            66.909090909091,
            29.440000000000,
            14.000000000000,
            14.000000000000,
            29.440000000000,
    };

    double q[] = { 0 - f[0], 1 - f[1]};

    struct ConvLFreeSet lfree;
    lfree.f = f;
    lfree.beta = beta;
    lfree.rays.nrays = 5;
    lfree.rays.values = rays;
    lfree.nrows = lfree.rays.dim = 2;

    double value;
    struct LP lp;

    rval = LP_open(&lp);
    abort_if(rval, "LP_open failed");

    rval = INFINITY_create_psi_lp(&lfree, &lp);
    abort_if(rval, "INFINITY_create_psi_lp failed");

    rval = INFINITY_psi(lfree.nrows, q, 1.0, &lp, &value);
    abort_if(rval, "GREDDY_ND_psi failed");
    EXPECT_NEAR(value, 1.0, 1e-6);

CLEANUP:
    if(rval) FAIL();
}

TEST(DISABLED_InfinityNDTest, generate_cut_test_1)
{
    int rval = 0;

    double r0[] = {  1.0,  1.0 };
    double r1[] = {  1.0, -1.0 };
    double r2[] = { -1.0, -1.0 };
    double r3[] = { -1.0,  1.0 };
    double r4[] = {  0.0,  1.0 };
    double r5[] = {  1.0,  0.0 };

    struct MultiRowModel model;
    CG_init_model(&model, 2, 6);
    LFREE_push_ray(&model.rays, r0);
    LFREE_push_ray(&model.rays, r1);
    LFREE_push_ray(&model.rays, r2);
    LFREE_push_ray(&model.rays, r3);
    LFREE_push_ray(&model.rays, r4);
    LFREE_push_ray(&model.rays, r5);
    model.f[0] = 0.5;
    model.f[1] = 0.5;

    struct ConvLFreeSet lfree;
    LFREE_init_conv(&lfree, 2, 6);

    rval = INFINITY_ND_generate_lfree(&model, &lfree);
    abort_if(rval, "INFINITY_ND_generate_lfree failed");

    EXPECT_NEAR(lfree.beta[0], 0.5, 1e-6);
    EXPECT_NEAR(lfree.beta[1], 0.5, 1e-6);
    EXPECT_NEAR(lfree.beta[2], 0.5, 1e-6);
    EXPECT_NEAR(lfree.beta[3], 0.5, 1e-6);
    EXPECT_NEAR(lfree.beta[4], 1.0, 1e-6);
    EXPECT_NEAR(lfree.beta[5], 1.0, 1e-6);

CLEANUP:
    if(rval) FAIL();
}

TEST(InfinityNDTest, generate_cut_test_2)
{
    int rval = 0;

    double r0[] = { 1.0,  0.0,  0.0 };
    double r1[] = {-1.0,  0.0,  0.0 };
    double r2[] = { 0.0,  1.0,  0.0 };
    double r3[] = { 0.0, -1.0,  0.0 };
    double r4[] = { 0.0,  0.0,  1.0 };
    double r5[] = { 0.0,  0.0, -1.0 };

    struct MultiRowModel model;
    CG_init_model(&model, 3, 6);
    LFREE_push_ray(&model.rays, r0);
    LFREE_push_ray(&model.rays, r1);
    LFREE_push_ray(&model.rays, r2);
    LFREE_push_ray(&model.rays, r3);
    LFREE_push_ray(&model.rays, r4);
    LFREE_push_ray(&model.rays, r5);
    model.f[0] = 0.75;
    model.f[1] = 0.75;
    model.f[2] = 0.75;

    struct ConvLFreeSet lfree;
    LFREE_init_conv(&lfree, 3, 6);

    rval = INFINITY_ND_generate_lfree(&model, &lfree);
    abort_if(rval, "INFINITY_ND_generate_lfree failed");

    EXPECT_NEAR(lfree.beta[0], 0.75, 1e-6);
    EXPECT_NEAR(lfree.beta[1], 2.25, 1e-6);
    EXPECT_NEAR(lfree.beta[2], 0.75, 1e-6);
    EXPECT_NEAR(lfree.beta[3], 2.25, 1e-6);
    EXPECT_NEAR(lfree.beta[4], 0.75, 1e-6);
    EXPECT_NEAR(lfree.beta[5], 2.25, 1e-6);

CLEANUP:
    CG_free_model(&model);
    if(rval) FAIL();
}

TEST(InfinityNDTest, scale_to_ahull_test)
{
    int rval = 0;

    double rays[] = 
    {
         0.0, 0.0, 1.0,
         0.0, 1.0, 0.0,
         1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0
    };

    int rx[] = { 1, 1, 1 , 0 };
    double beta[] = { INFINITY, INFINITY, INFINITY, INFINITY };
    double epsilon = 1.0;

    double d1[] = {  1.0,  1.0,  1.0 };
    double d2[] = {  2.0,  2.0,  0.0 };
    double d3[] = { -1.0, -1.0, -1.0 };

    double alpha;

    rval = scale_to_ahull(3, 4, rays, rx, beta, epsilon, d1, &alpha);
    abort_if(rval, "scale_to_ahull failed");
    EXPECT_DOUBLE_EQ(1 / 3.0, alpha);

    rval = scale_to_ahull(3, 4, rays, rx, beta, epsilon, d2, &alpha);
    abort_if(rval, "scale_to_ahull failed");
    EXPECT_DOUBLE_EQ(0.25, alpha);

    rval = scale_to_ahull(3, 4, rays, rx, beta, epsilon, d3, &alpha);
    abort_if(rval, "scale_to_ahull failed");
    EXPECT_EQ(INFINITY, alpha);

CLEANUP:
    if(rval) FAIL();
}
