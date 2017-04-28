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

extern "C" {
#include <multirow/util.h>
#include <multirow/double.h>
#include <multirow/lfree2d.h>
#include <lifting/lifting.h>
#include <lifting/lifting-mip.h>
#include <math.h>
}

#define E 1e-6

TEST(Lifting2DTest, lift_fixed_test_2)
{
    int rval = 0;

    double opt;

    int n_halfspaces = 3;
    double halfspaces[] = {
            -1.12796209,  -2.25592417,
            2.38125329,  11.19606810,
            7.35264174,   6.22467965,
    };

    double r[] = { 0.48231629, 0.70551355 };

    rval = LIFTING_2D_lift_fixed(n_halfspaces, halfspaces, r, -1, &opt);
    abort_if(rval, "LIFTING_lift_fixed failed");
    EXPECT_NEAR(opt, 1.24826670, E);

CLEANUP:
    if(rval) FAIL();
}

TEST(Lifting2DTest, psi_test)
{
    int rval = 0;

    int n_halfspaces = 5;
    double halfspaces[] = {
        -1 / 2.0, 0,
        -1 / 2.0, -1 / 2.0,
        -1 / 4.0, 1 / 2.0,
        -1 / 6.0, -5 / 6.0,
         5 / 16.0, 1 / 8.0
    };

    double r1[] = { 0, 0 };
    double r2[] = { 2, 3 };
    double r3[] = { 1 / 3.0, 1 / 2.0 };
    double r4[] = { -4.0, 0 };
    double r5[] = { 1, -1/2.0 };
    double r6[] = { -1/2.0, 1/2.0 };
    double r7[] = { 1/2.0, 1/2.0 };
    double value;

    rval = LIFTING_2D_psi(n_halfspaces, halfspaces, r1, &value);
    abort_if(rval, "LIFTING_2D_psi failed");
    EXPECT_NEAR(value, 0.0, E);

    rval = LIFTING_2D_psi(n_halfspaces, halfspaces, r2, &value);
    abort_if(rval, "LIFTING_2D_psi failed");
    EXPECT_NEAR(value, 1.0, E);

    rval = LIFTING_2D_psi(n_halfspaces, halfspaces, r3, &value);
    abort_if(rval, "LIFTING_2D_psi failed");
    EXPECT_NEAR(value, 1 / 6.0, E);

    rval = LIFTING_2D_psi(n_halfspaces, halfspaces, r4, &value);
    abort_if(rval, "LIFTING_2D_psi failed");
    EXPECT_NEAR(value, 2.0, E);

    rval = LIFTING_2D_psi(n_halfspaces, halfspaces, r5, &value);
    abort_if(rval, "LIFTING_2D_psi failed");
    EXPECT_NEAR(value, 1 / 4.0, E);

    rval = LIFTING_2D_psi(n_halfspaces, halfspaces, r6, &value);
    abort_if(rval, "LIFTING_2D_psi failed");
    EXPECT_NEAR(value, 0.375, E);

    rval = LIFTING_2D_psi(n_halfspaces, halfspaces, r7, &value);
    abort_if(rval, "LIFTING_2D_psi failed");
    EXPECT_NEAR(value, 7 / 32.0, E);

CLEANUP:
    if(rval) FAIL();
}

TEST(Lifting2DTest, psi_test_2)
{
    int rval = 0;

    LFreeSet2D set;

    set.f[0] = 0;
    set.f[1] = 0.359;

    double vertices[] = {
        -1.700870, -0.796700,
         0.068700,  0.032200,
         1.106900,  1.111100,
        -0.131990,  0.986800,
    };

    double lattice_points[] = {
         0.000000, 0.000000,
         1.000000, 1.000000,
         0.000000, 1.000000,
        -1.000000, 0.000000,
    };


    int ub[] = { 10, 10 };
    int lb[] = { -12, -12 };

    double r0[] = { 0.37, 0.19 };
    double value;

    rval = LFREE_2D_init(&set, 100, 100, 100);
    abort_if(rval, "LFREE_2D_init failed");

    set.n_vertices = 4;
    set.vertices = vertices;

    set.n_lattice_points = 4;
    set.lattice_points = lattice_points;

    rval = LFREE_2D_compute_halfspaces(&set);
    abort_if(rval, "LFREE_2D_compute_halfspaces failed");

    rval = LFREE_2D_print_set(&set);
    abort_if(rval, "LFREE_2D_print_set failed");

    rval = LIFTING_2D_naive(set.n_halfspaces, set.halfspaces, r0, lb, ub, &value);
    abort_if(rval, "LIFTING_2D_psi failed");

    EXPECT_NEAR(value, 0.48846868, E);

CLEANUP:
    if(rval) FAIL();
}

TEST(Lifting2DTest, naive_test_1)
{
    int rval = 0;

    int n_halfspaces = 5;
    double halfspaces[] = {
        -1 / 2.0, 0,
        -1 / 2.0, -1 / 2.0,
        -1 / 4.0, 1 / 2.0,
        -1 / 6.0, -5 / 6.0,
         5 / 16.0, 1 / 8.0
    };

    double r1[] = { 1/2.0, 1/2.0 };
    double r2[] = { 0, 1/2.0 };
    double r3[] = { 1/3.0, 2/3.0 };
    double r4[] = { 2/5.0, 3/7.0 };
    int lb[] = { -10, -10 };
    int ub[] = {  10,  10 };

    double value;

    rval = LIFTING_2D_naive(n_halfspaces, halfspaces, r1, lb, ub, &value);
    abort_if(rval, "LIFTING_2D_naive failed");
    EXPECT_NEAR(value, 0.21875, E);

    rval = LIFTING_2D_naive(n_halfspaces, halfspaces, r2, lb, ub, &value);
    abort_if(rval, "LIFTING_2D_naive failed");
    EXPECT_NEAR(value, 0.25, E);

    rval = LIFTING_2D_naive(n_halfspaces, halfspaces, r3, lb, ub, &value);
    abort_if(rval, "LIFTING_2D_naive failed");
    EXPECT_NEAR(value, 0.222222, E);

    rval = LIFTING_2D_naive(n_halfspaces, halfspaces, r4, lb, ub, &value);
    abort_if(rval, "LIFTING_2D_naive failed");
    EXPECT_NEAR(value, 0.178571, E);

CLEANUP:
    if(rval) FAIL();
}

TEST(Lifting2DTest, naive_test_2)
{
    int rval = 0;

    int n_halfspaces = 3;
    double halfspaces[] = {
            -1.12796209,  -2.25592417,
            2.38125329,  11.19606810,
            7.35264174,   6.22467965,
    };

    double r[] = { 0.48231629, 0.70551355 };
    int lb[] = { -50, -50 };
    int ub[] = {  50,  50 };

    double value;

    rval = LIFTING_2D_naive(n_halfspaces, halfspaces, r, lb, ub, &value);
    abort_if(rval, "LIFTING_2D_naive failed");
    EXPECT_NEAR(value, 1.24826671, E);

CLEANUP:
    if(rval) FAIL();
}

TEST(Lifting2DTest, bound_test_1)
{
    int rval = 0;

    int n_halfspaces = 5;
    double halfspaces[] = {
            -1 / 2.0, 0,
            -1 / 2.0, -1 / 2.0,
            -1 / 4.0, 1 / 2.0,
            -1 / 6.0, -5 / 6.0,
            5 / 16.0, 1 / 8.0
    };

    double r1[] = { 1/2.0, 1/2.0 };
    double r2[] = { 0, 1/2.0 };
    double r3[] = { 1/3.0, 2/3.0 };
    double r4[] = { 2/5.0, 3/7.0 };

    double value;

    rval = LIFTING_2D_bound(n_halfspaces, halfspaces, r1, &value);
    abort_if(rval, "LIFTING_2D_bound failed");
    EXPECT_NEAR(value, 0.21875, E);

    rval = LIFTING_2D_bound(n_halfspaces, halfspaces, r2, &value);
    abort_if(rval, "LIFTING_2D_bound failed");
    EXPECT_NEAR(value, 0.25, E);

    rval = LIFTING_2D_bound(n_halfspaces, halfspaces, r3, &value);
    abort_if(rval, "LIFTING_2D_bound failed");
    EXPECT_NEAR(value, 0.222222, E);

    rval = LIFTING_2D_bound(n_halfspaces, halfspaces, r4, &value);
    abort_if(rval, "LIFTING_2D_bound failed");
    EXPECT_NEAR(value, 0.178571, E);

CLEANUP:
    if(rval) FAIL();
}

TEST(Lifting2DTest, bound_test_2)
{
    int rval = 0;

    int n_halfspaces = 3;
    double halfspaces[] = {
        -1.12796209,  -2.25592417,
         2.38125329,  11.19606810,
         7.35264174,   6.22467965,
    };

    double r[] = { 0.48231629, 0.70551355 };

    double value;

    rval = LIFTING_2D_bound(n_halfspaces, halfspaces, r, &value);
    abort_if(rval, "LIFTING_2D_bound failed");
    EXPECT_NEAR(value, 1.24826671, E);

CLEANUP:
    if(rval) FAIL();
}

TEST(Lifting2DTest, mip_test_1)
{
    int rval = 0;

    int n_halfspaces = 5;
    double halfspaces[] = {
            -1 / 2.0, 0,
            -1 / 2.0, -1 / 2.0,
            -1 / 4.0, 1 / 2.0,
            -1 / 6.0, -5 / 6.0,
            5 / 16.0, 1 / 8.0
    };

    double r1[] = { 1/2.0, 1/2.0 };
    double r2[] = { 0, 1/2.0 };
    double r3[] = { 1/3.0, 2/3.0 };
    double r4[] = { 2/5.0, 3/7.0 };

    double value;

    rval = LIFTING_2D_mip(n_halfspaces, halfspaces, r1, &value);
    abort_if(rval, "LIFTING_2D_mip failed");
    EXPECT_NEAR(value, 0.21875, E);

    rval = LIFTING_2D_mip(n_halfspaces, halfspaces, r2, &value);
    abort_if(rval, "LIFTING_2D_mip failed");
    EXPECT_NEAR(value, 0.25, E);

    rval = LIFTING_2D_mip(n_halfspaces, halfspaces, r3, &value);
    abort_if(rval, "LIFTING_2D_mip failed");
    EXPECT_NEAR(value, 0.222222, E);

    rval = LIFTING_2D_mip(n_halfspaces, halfspaces, r4, &value);
    abort_if(rval, "LIFTING_2D_mip failed");
    EXPECT_NEAR(value, 0.178571, E);

    CLEANUP:
    if(rval) FAIL();
}

TEST(Lifting2DTest, verify_test_1)
{
    int rval = 0;

    double vertices[] = {
         0.361588, -0.411851,
         2.550110,  1.000000,
         1.009740,  1.000000,
        -0.836467,  0.952742,
    };

    double lattice_points[] = {
         1.000000,  0.000000,
         2.000000,  1.000000,
         1.000000,  1.000000,
         0.000000,  0.000000,
    };

    LFreeSet2D set;

    rval = LFREE_2D_init(&set, 4, 4, 4);
    abort_if(rval, "LFREE_2D_init failed");

    set.n_vertices = 4;
    set.vertices = vertices;

    set.n_lattice_points = 4;
    set.lattice_points = lattice_points;

    set.f[0] = 0.918857;
    set.f[1] = 0.952742;

    rval = LFREE_2D_compute_halfspaces(&set);
    abort_if(rval, "LFREE_2D_compute_halfspaces failed");

    rval = LIFTING_2D_verify(&set);
    EXPECT_EQ(rval, 1);
    rval = 0;

CLEANUP:
    if(rval) FAIL();
}

TEST(LFreeSetTest, read_next_test)
{
    int rval = 0;

    FILE *triangles = fopen("../lifting/library/tests/fixtures/triangles.txt", "r");
    abort_if(!triangles, "could not read triangles.txt");

    LFreeSet2D set;
    LFREE_2D_init(&set, 100, 100, 100);

    rval = LFREE_2D_read_next(triangles, &set);
    abort_if(rval, "LFREE_2D_read_next failed");

    EXPECT_NEAR(set.f[0], 0.7875, E);
    EXPECT_NEAR(set.f[1], 1.7875, E);

    EXPECT_NEAR(set.vertices[0], 1, E);
    EXPECT_NEAR(set.vertices[1], 2, E);

    EXPECT_NEAR(set.vertices[2], 0.4875, E);
    EXPECT_NEAR(set.vertices[3], 2.0875, E);

    EXPECT_NEAR(set.vertices[4], -0.0625, E);
    EXPECT_NEAR(set.vertices[5], -16.0625, E);

    EXPECT_NEAR(set.lattice_points[0], 1, E);
    EXPECT_NEAR(set.lattice_points[1], 2, E);

    EXPECT_NEAR(set.lattice_points[2], 0, E);
    EXPECT_NEAR(set.lattice_points[3], -14, E);

    EXPECT_NEAR(set.lattice_points[4], 0, E);
    EXPECT_NEAR(set.lattice_points[5], -15, E);

    rval = LFREE_2D_compute_halfspaces(&set);
    abort_if(rval, "LFREE_2D_compute_halfspaces failed");

    EXPECT_NEAR(set.halfspaces[0], 0.686274, E);
    EXPECT_NEAR(set.halfspaces[1], 4.019607, E);

    EXPECT_NEAR(set.halfspaces[2], -3.235294, E);
    EXPECT_NEAR(set.halfspaces[3],  0.098039, E);

    EXPECT_NEAR(set.halfspaces[4],  5.0, E);
    EXPECT_NEAR(set.halfspaces[5], -0.294117, E);

CLEANUP:
    if(rval) FAIL();
}

TEST(LFreeSetTest, read_next_quadrilateral_test)
{
    int rval = 0;

    FILE *file = fopen("../lifting/library/tests/fixtures/quads.txt", "r");
    abort_if(!file, "could not read quads.txt");

    LFreeSet2D set;
    LFREE_2D_init(&set, 100, 100, 100);

    rval = LFREE_2D_read_next(file, &set);
    abort_if(rval, "LFREE_2D_read_next failed");

    EXPECT_NEAR(set.f[0], 0.5, E);
    EXPECT_NEAR(set.f[1], 0.5, E);

    EXPECT_NEAR(set.vertices[0], -0.5, E);
    EXPECT_NEAR(set.vertices[1], 0.5, E);

    EXPECT_NEAR(set.vertices[2], 0.5, E);
    EXPECT_NEAR(set.vertices[3], 1.5, E);

    EXPECT_NEAR(set.vertices[4], 1.5, E);
    EXPECT_NEAR(set.vertices[5], 0.5, E);

    EXPECT_NEAR(set.vertices[6], 0.5, E);
    EXPECT_NEAR(set.vertices[7], -0.5, E);

    EXPECT_NEAR(set.lattice_points[0], 0, E);
    EXPECT_NEAR(set.lattice_points[1], 0, E);

    EXPECT_NEAR(set.lattice_points[2], 0, E);
    EXPECT_NEAR(set.lattice_points[3], 1, E);

    EXPECT_NEAR(set.lattice_points[4], 1, E);
    EXPECT_NEAR(set.lattice_points[5], 1, E);

    EXPECT_NEAR(set.lattice_points[6], 1, E);
    EXPECT_NEAR(set.lattice_points[7], 0, E);

    rval = LFREE_2D_compute_halfspaces(&set);
    abort_if(rval, "LFREE_2D_compute_halfspaces failed");

    EXPECT_NEAR(set.halfspaces[0], -1, E);
    EXPECT_NEAR(set.halfspaces[1], 1, E);

    EXPECT_NEAR(set.halfspaces[2], 1, E);
    EXPECT_NEAR(set.halfspaces[3], 1, E);

    EXPECT_NEAR(set.halfspaces[4], 1, E);
    EXPECT_NEAR(set.halfspaces[5], -1, E);

    EXPECT_NEAR(set.halfspaces[6], -1, E);
    EXPECT_NEAR(set.halfspaces[7], -1, E);

CLEANUP:
    if(rval) FAIL();
}

//TEST(LFreeSetTest, get_bounding_box_test)
//{
//    int rval = 0;
//
//    FILE *triangles = fopen("../tests/fixtures/triangles.txt", "r");
//    abort_if(!triangles, "could not read triangles.txt");
//
//    LFreeSet2D set;
//    LFREE_2D_init(&set, 100, 100, 100);
//
//    int ub[2], lb[2];
//
//    rval = LFREE_2D_read_next(triangles, &set);
//    abort_if(rval, "LFREE_2D_read_next failed");
//
//    rval = LFREE_2D_get_bounding_box(&set, lb, ub);
//    abort_if(rval, "LFREE_2D_get_bounding_box failed");
//
//    EXPECT_NEAR(lb[0], -1.0, E);
//    EXPECT_NEAR(ub[0], 1.0, E);
//
//    EXPECT_NEAR(lb[1], -18.0, E);
//    EXPECT_NEAR(ub[1], 2, E);
//
//CLEANUP:
//    if(rval) FAIL();
//}

TEST(LFreeSetTest, preprocess_test)
{
    int rval = 0;

    double vertices[] = {
            100/11.0 ,167/33.0,
             10/11.0 , 97/33.0,
             20/11.0 , 39/11.0
    };

    double lattice_points[] = {
            1, 3,
            4, 4,
            5, 4
    };

    double pre_m[4];
    double center;

    struct LFreeSet2D set;
    set.f[0] = 9 / 2.0;
    set.f[1] = 4.0;

    set.n_vertices = 3;
    set.vertices = vertices;

    set.n_lattice_points = 3;
    set.lattice_points = lattice_points;

    rval = LFREE_2D_preprocess(&set, pre_m, &center);
    abort_if(rval, "LFREE_2D_preprocess failed");

    EXPECT_NEAR(set.f[0], 1 / 2.0, E);
    EXPECT_NEAR(set.f[1], 1 / 2.0, E);

    EXPECT_NEAR(set.vertices[0], 21 / 11.0, E);
    EXPECT_NEAR(set.vertices[1], 5 / 33.0, E);

    EXPECT_NEAR(set.vertices[2],   1 / 11.0, E);
    EXPECT_NEAR(set.vertices[3], -5 / 33.0, E);

    EXPECT_NEAR(set.vertices[4], -9 / 11.0, E);
    EXPECT_NEAR(set.vertices[5], 15 / 11.0, E);

    EXPECT_NEAR(center, 20 / 33.0, E);

CLEANUP:
    if(rval) FAIL();

}TEST(LFreeSetTest, preprocess_test_3)
{
    int rval = 0;

    double vertices[] = {
             1.0, -4.0,
             1.0,  2.0,
            -0.2,  0.8
    };

    double lattice_points[] = {
            1, 0,
            0, 1,
            0, 0
    };

    double pre_m[4];
    double center;

    struct LFreeSet2D set;
    set.f[0] = 1 / 2.0;
    set.f[1] = 1 / 2.0;

    set.n_vertices = 3;
    set.vertices = vertices;

    set.n_lattice_points = 3;
    set.lattice_points = lattice_points;

    rval = LFREE_2D_preprocess(&set, pre_m, &center);
    abort_if(rval, "LFREE_2D_preprocess failed");

    EXPECT_NEAR(set.vertices[0], -4, E);
    EXPECT_NEAR(set.vertices[1], 1, E);

    EXPECT_NEAR(set.vertices[2], 2, E);
    EXPECT_NEAR(set.vertices[3], 1, E);

    EXPECT_NEAR(set.vertices[4], 0.8, E);
    EXPECT_NEAR(set.vertices[5], -0.2, E);

    EXPECT_NEAR(center, 0.4, E);

CLEANUP:
    if(rval) FAIL();
}

TEST(LFreeSetTest, preprocess_test_2)
{
    int rval = 0;

    double vertices[] = {
        22, 69/7.0,
        -3, -11/7.0,
        -8, -26/7.0
    };

    double lattice_points[] = {
        -4, -2,
        -2, -1,
        7, 3
    };

    double r[] = { 2/3.0, 2/6.0 };
    double r0, r1;
    double value;

    double pre_m[4];

    struct LFreeSet2D set;
    set.f[0] = 2 / 3.0;
    set.f[1] = 1 / 6.0;

    rval = LFREE_2D_init(&set, 10, 10, 10);
    abort_if(rval, "LFREE_2D_init failed");

    set.n_vertices = 3;
    set.vertices = vertices;

    set.n_lattice_points = 3;
    set.lattice_points = lattice_points;

    rval = LFREE_2D_print_set(&set);
    abort_if(rval, "LFREE_2D_print_set failed");

    rval = LFREE_2D_compute_halfspaces(&set);
    abort_if(rval, "LFREE_2D_compute_halfspaces failed");

    rval = LIFTING_2D_bound(set.n_halfspaces, set.halfspaces, r, &value);
    abort_if(rval, "LIFTING_2D_bound failed");

    log_debug("%.6lf\n", value);

CLEANUP:
    if(rval) FAIL();
}


TEST(LFreeSetTest, preprocess_test_4)
{
    int rval = 0;

    double vertices[] = {
            -5/2.0, -1,
            17/4.0, 2,
            7/2.0, 2,
            -13/4.0, -1
    };

    double lattice_points[] = {
            2, 1,
            -3, -1,
            -1, 0,
            4, 2
    };

    double pre_m[4];
    double center;

    struct LFreeSet2D set;
    set.f[0] = 1 / 2.0;
    set.f[1] = 1 / 2.0;

    set.n_vertices = 4;
    set.vertices = vertices;

    set.n_lattice_points = 4;
    set.lattice_points = lattice_points;

    rval = LFREE_2D_preprocess(&set, pre_m, &center);
    abort_if(rval, "LFREE_2D_preprocess failed");

    LFREE_2D_print_set(&set);

    EXPECT_NEAR(set.vertices[0], -1, E);
    EXPECT_NEAR(set.vertices[1], 1/2.0, E);

    EXPECT_NEAR(set.vertices[2], 1/2.0, E);
    EXPECT_NEAR(set.vertices[3], -1/4.0, E);

    EXPECT_NEAR(set.vertices[4], 2, E);
    EXPECT_NEAR(set.vertices[5], 1/2.0, E);

    EXPECT_NEAR(set.vertices[6], 1/2.0, E);
    EXPECT_NEAR(set.vertices[7], 5/4.0, E);

    EXPECT_NEAR(pre_m[0], -2, E);
    EXPECT_NEAR(pre_m[1], -1, E);
    EXPECT_NEAR(pre_m[2], 5, E);
    EXPECT_NEAR(pre_m[3], 2, E);

    EXPECT_NEAR(center, 0.5, E);

    CLEANUP:
    if(rval) FAIL();
}
