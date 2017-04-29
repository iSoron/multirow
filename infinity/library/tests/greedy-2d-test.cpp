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
using namespace std;

#define TEST_SOURCE

extern "C" {
#include <multirow/util.h>
#include <infinity/greedy-2d.h>
#include "../src/greedy-2d.c"
}

#define BOUNDS_EPSILON 0.01

TEST(Greedy2DTest, test_generate_cut_1)
{
    int rval = 0;
    double bounds[100];
    double f[] = {1 / 4.0, 3 / 4.0};
    double rays[] = {-2 / 5.0, 5 / 7.0, 0.0, 1.0, 1.0, 1.0, 4 / 5.0, -2 / 3.0,
            -1.0, 0.0};

    rval = GREEDY_2D_generate_cut(rays, 5, f, bounds);
    abort_if(rval, "GREEDY_2D_generate_cut failed");

    EXPECT_NEAR(23 / 50.0, bounds[0], BOUNDS_EPSILON);
    EXPECT_NEAR(23 / 42.0, bounds[1], BOUNDS_EPSILON);
    EXPECT_NEAR(9 / 11.0, bounds[2], BOUNDS_EPSILON);
    EXPECT_NEAR(9 / 11.0, bounds[3], BOUNDS_EPSILON);
    EXPECT_NEAR(23 / 50.0, bounds[4], BOUNDS_EPSILON);

    CLEANUP:
    if (rval) FAIL();
}

TEST(Greedy2DTest, test_generate_cut_2)
{
    int rval = 0;
    double bounds[100];
    double f[] = {1 / 2.0, 1 / 2.0};
    double rays[] = {
            -1.0, -1.0,
            -1.0, 1.0,
            1.0, 1.0,
            1.0, 0.0,
            1.0, -1.0
    };

    rval = GREEDY_2D_generate_cut(rays, 5, f, bounds);
    abort_if(rval, "GREEDY_2D_generate_cut failed");

    EXPECT_NEAR(0.5, bounds[0], BOUNDS_EPSILON);
    EXPECT_NEAR(0.5, bounds[1], BOUNDS_EPSILON);
    EXPECT_NEAR(0.5, bounds[2], BOUNDS_EPSILON);
    EXPECT_EQ(GREEDY_BIG_E, bounds[3]);
    EXPECT_NEAR(0.5, bounds[4], BOUNDS_EPSILON);

    CLEANUP:
    if (rval) FAIL();
}

TEST(Greedy2DTest, test_generate_cut_3)
{
    int rval = 0;
    double bounds[100];
    double f[] = {5 / 22.0, 0.0};
    double rays[] = {-1 / 22.0, 0.0, 0.0, 1 / 18.0, 1 / 22.0, 0.0};

    rval = GREEDY_2D_generate_cut(rays, 3, f, bounds);
    abort_if(rval, "GREEDY_2D_generate_cut failed");

    EXPECT_NEAR(5.0, bounds[0], BOUNDS_EPSILON);
    EXPECT_NEAR(17.0, bounds[2], BOUNDS_EPSILON);
    EXPECT_EQ(GREEDY_BIG_E, bounds[1]);

    CLEANUP:
    if (rval) FAIL();
}

TEST(Greedy2DTest, scale_to_chull_test)
{
    int rval = 0;

    double rays[] = {
            0, 1,
            1, 1,
            0.5, 0.25,
            1, 0,
            -1, -1,
            -0.25, 0
    };

    double scale[6];
    int nrays = 6;

    rval = scale_to_chull(rays, nrays, scale);
    abort_if(rval, "scale_to_chull failed");

    EXPECT_NEAR(rays[0], 0.0, BOUNDS_EPSILON);
    EXPECT_NEAR(rays[1], 1.0, BOUNDS_EPSILON);

    EXPECT_NEAR(rays[2], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(rays[3], 1.0, BOUNDS_EPSILON);

    EXPECT_NEAR(rays[4], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(rays[5], 0.5, BOUNDS_EPSILON);

    EXPECT_NEAR(rays[6], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(rays[7], 0.0, BOUNDS_EPSILON);

    EXPECT_NEAR(rays[8],-1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(rays[9],-1.0, BOUNDS_EPSILON);

    EXPECT_NEAR(rays[10],-0.5, BOUNDS_EPSILON);
    EXPECT_NEAR(rays[11], 0.0, BOUNDS_EPSILON);

    EXPECT_NEAR(scale[0], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(scale[1], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(scale[2], 2.0, BOUNDS_EPSILON);
    EXPECT_NEAR(scale[3], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(scale[4], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(scale[5], 2.0, BOUNDS_EPSILON);

CLEANUP:
    if(rval) FAIL();
}

TEST(Greedy2DTest, scale_to_chull_test_2)
{
    int rval = 0;

    double rays[] = {
            1, 1,
            0.5, 0.25,
            1, 0,
    };

    double scale[3];
    int nrays = 3;

    rval = scale_to_chull(rays, nrays, scale);
    abort_if(rval, "scale_to_chull failed");

    EXPECT_NEAR(rays[0], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(rays[1], 1.0, BOUNDS_EPSILON);

    EXPECT_NEAR(rays[2], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(rays[3], 0.5, BOUNDS_EPSILON);

    EXPECT_NEAR(rays[4], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(rays[5], 0.0, BOUNDS_EPSILON);

    EXPECT_NEAR(scale[0], 1.0, BOUNDS_EPSILON);
    EXPECT_NEAR(scale[1], 2.0, BOUNDS_EPSILON);
    EXPECT_NEAR(scale[2], 1.0, BOUNDS_EPSILON);

    CLEANUP:
    if(rval) FAIL();
}

TEST(Greedy2DTest, find_containing_cone_test)
{
    int rval = 0;

    double rays[] = {-1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, -1.0, -1.0, -1.0};
    double p1[] = {-1.0, 1.0};
    double p2[] = {0.0, 1.0};
    double p3[] = {1.0, 0.0};
    double p4[] = {-1.0, 0.0};
    int index1, index2;
    double lambda1, lambda2;

    rval = find_containing_cone(rays, 5, p1, &index1, &index2, &lambda1, &lambda2);
    abort_if(rval, "find_containing_cone failed");

    EXPECT_EQ(0, index1);
    EXPECT_EQ(1, index2);

    rval = find_containing_cone(rays, 5, p2, &index1, &index2, &lambda1, &lambda2);
    abort_if(rval, "find_containing_cone failed");

    EXPECT_EQ(0, index1);
    EXPECT_EQ(1, index2);

    rval = find_containing_cone(rays, 5, p3, &index1, &index2, &lambda1, &lambda2);
    abort_if(rval, "find_containing_cone failed");

    EXPECT_EQ(1, index1);
    EXPECT_EQ(2, index2);

    rval = find_containing_cone(rays, 5, p4, &index1, &index2, &lambda1, &lambda2);
    abort_if(rval, "find_containing_cone failed");

    EXPECT_EQ(4, index1);
    EXPECT_EQ(0, index2);

    CLEANUP:
    if(rval) FAIL();
}

TEST(Greedy2DTest, find_containing_cone_test_2)
{
    int rval = 0;

    double rays[] = {0.0, 0.5, 1.0, 0.0, -0.5, 0.5};
    double p[] = {0.5, -0.5};
    int index1, index2;
    double lambda1, lambda2;

    rval = find_containing_cone(rays, 3, p, &index1, &index2, &lambda1, &lambda2);
    abort_if(rval, "find_containing_cone failed");

    EXPECT_EQ(-1, index1);
    EXPECT_EQ(-1, index2);

    CLEANUP:
    if(rval) FAIL();
}

TEST(Greedy2DTest, find_containing_cone_test_3)
{
    int rval = 0;

    double rays[] = {
            -1, -0.1,
             1, -0.1,
             1, -0.5,
            -1, -0.5,
    };

    double p[] = { 0, -2 };
    int index1, index2;
    double lambda1, lambda2;

    rval = find_containing_cone(rays, 4, p, &index1, &index2, &lambda1, &lambda2);
    abort_if(rval, "find_containing_cone failed");

    EXPECT_EQ(2, index1);
    EXPECT_EQ(3, index2);

 CLEANUP:
    if(rval) FAIL();
}

//TEST(Greedy2DTest, test_generate_cut_4)
//{
//    int rval = 0;
//    double bounds[100];
//    double f[] = {5 / 22.0, 10 / 19.0};
//    double rays[] = {0, -1 / 38.0, -1 / 22.0, -1 / 38.0, 0, 1 / 38.0, -1 / 22.0,
//            0, 1 / 22.0, 0, 1 / 22.0, 1 / 38.0};
//
//    rval = GREEDY_2D_generate_cut(rays, 6, f, bounds);
//    abort_if(rval, "GREEDY_2D_generate_cut failed");
//
//    EXPECT_NEAR(20.0, bounds[0], BOUNDS_EPSILON);
//    EXPECT_NEAR(20.0, bounds[1], BOUNDS_EPSILON);
//    EXPECT_NEAR(18.0, bounds[2], BOUNDS_EPSILON);
//    EXPECT_NEAR(18.0, bounds[5], BOUNDS_EPSILON);
//    EXPECT_EQ(GREEDY_BIG_E, bounds[3]);
//    EXPECT_EQ(GREEDY_BIG_E, bounds[4]);
//
//    CLEANUP:
//    if (rval) FAIL();
//}
//
//TEST(Greedy2DTest, test_generate_cut_5)
//{
//    int rval = 0;
//    double bounds[100];
//    double f[] = {0.22727272727272729291, 0.52631578947368418131};
//    double rays[] = {0.00000000000000000000, -0.02631578947368420907,
//            -0.04545454545454545581, -0.02631578947368420907,
//            0.00000000000000000000, 0.02631578947368420907,
//            -0.04545454545454545581, 0.00000000000000000000,
//            0.04545454545454545581, 0.00000000000000000000,
//            0.04545454545454545581, 0.02631578947368420907};
//
//    rval = GREEDY_2D_generate_cut(rays, 6, f, bounds);
//    abort_if(rval, "GREEDY_2D_generate_cut failed");
//
//    EXPECT_NEAR(20.0, bounds[0], BOUNDS_EPSILON);
//    EXPECT_NEAR(20.0, bounds[1], BOUNDS_EPSILON);
//    EXPECT_NEAR(18.0, bounds[2], BOUNDS_EPSILON);
//    EXPECT_NEAR(18.0, bounds[5], BOUNDS_EPSILON);
//    EXPECT_EQ(GREEDY_BIG_E, bounds[3]);
//    EXPECT_EQ(GREEDY_BIG_E, bounds[4]);
//
//    CLEANUP:
//    if (rval) FAIL();
//}
//
//
//
//TEST(Greedy2DTest, get_peak_ray_test_1)
//{
//    int rval = 0;
//
//    double rays[] = {-2.0, 0.0, -1.0, 1.0, 1.0, 1.0, 2.0, 0.0, 1.0, -1.0, -1.0, -1.0 };
//    int nrays = 6;
//
//    double normal1[] = { 0.0,  1.0 };
//    double normal2[] = { 1.0,  1.0 };
//    double normal3[] = { 0.0, -1.0 };
//    double normal4[] = {-1.0, -1.0 };
//    int index;
//
//    rval = get_peak_ray(rays, nrays, normal1, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_TRUE(index == 1 || index == 2);
//
//    rval = get_peak_ray(rays, nrays, normal3, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_TRUE(index == 4 || index == 5);
//
//    rval = get_peak_ray(rays, nrays, normal2, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_EQ(2, index);
//
//    rval = get_peak_ray(rays, nrays, normal4, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_EQ(5, index);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(Greedy2DTest, get_peak_ray_test_2)
//{
//    int rval = 0;
//
//    double rays[] = { -1e100, 0, 0, 1, 1e100, 0, 1, -1, -1,-1 };
//    int nrays = 5;
//
//    double normal1[] = { 0.0, 1.0 };
//    double normal2[] = { 1.0, 1.0 };
//    double normal3[] = { 1.0, 0.0 };
//    double normal4[] = {  0.0, -1.0 };
//    double normal5[] = { -1.0, -1.0 };
//    double normal6[] = { -1.0,  0.0 };
//    int index;
//
//    rval = get_peak_ray(rays, nrays, normal1, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_EQ(1, index);
//
//    rval = get_peak_ray(rays, nrays, normal4, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_EQ(3, index);
//
//    rval = get_peak_ray(rays, nrays, normal2, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_EQ(2, index);
//
//    rval = get_peak_ray(rays, nrays, normal5, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_EQ(0, index);
//
//    rval = get_peak_ray(rays, nrays, normal3, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_EQ(2, index);
//
//    rval = get_peak_ray(rays, nrays, normal6, &index);
//    abort_if(rval, "get_peak_ray failed");
//    EXPECT_EQ(0, index);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(Greedy2DTest, get_peak_ray_test_3)
//{
//    int rval = 0;
//
//    double rays[] = { -1, 1, 0, 1, 1, 1};
//    int nrays = 3;
//
//    double normal1[] = { 0.0, -1.0 };
//    int index;
//
//    rval = get_peak_ray(rays, nrays, normal1, &index);
//    EXPECT_EQ(ERR_NOT_FOUND, rval);
//    rval = 0;
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(Greedy2DTest, nearest_lattice_point_test_1)
//{
//    int rval = 0;
//
//    double z0a[] = {100, 1};
//    double z0b[] = {-100, 1};
//    double g[] = {-1, 0};
//
//    double x1[] = {-1.34, 1};
//    double x2[] = {-3, 1};
//
//    double z1[2];
//    double z2[2];
//
//    rval = nearest_lattice_points(g, z0a, x1, z1, z2);
//    abort_if(rval, "nearest_lattice_points failed");
//    EXPECT_DOUBLE_EQ(-1, z1[0]);
//    EXPECT_DOUBLE_EQ(1, z1[1]);
//    EXPECT_DOUBLE_EQ(-2, z2[0]);
//    EXPECT_DOUBLE_EQ(1, z2[1]);
//
//    rval = nearest_lattice_points(g, z0b, x1, z1, z2);
//    abort_if(rval, "nearest_lattice_points failed");
//    EXPECT_DOUBLE_EQ(-1, z1[0]);
//    EXPECT_DOUBLE_EQ(1, z1[1]);
//    EXPECT_DOUBLE_EQ(-2, z2[0]);
//    EXPECT_DOUBLE_EQ(1, z2[1]);
//
//    rval = nearest_lattice_points(g, z0a, x2, z1, z2);
//    abort_if(rval, "nearest_lattice_points failed");
//    EXPECT_DOUBLE_EQ(-3, z1[0]);
//    EXPECT_DOUBLE_EQ(1, z1[1]);
//    EXPECT_DOUBLE_EQ(-3, z2[0]);
//    EXPECT_DOUBLE_EQ(1, z2[1]);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(Greedy2DTest, nearest_lattice_point_test_2)
//{
//    int rval = 0;
//
//    double z0[] = {1, 1};
//    double g[] = {1, 2};
//    double x[] = {3.35, 5.7};
//
//    double z1[2];
//    double z2[2];
//
//    rval = nearest_lattice_points(g, z0, x, z1, z2);
//    abort_if(rval, "nearest_lattice_points failed");
//    EXPECT_DOUBLE_EQ(3, z1[0]);
//    EXPECT_DOUBLE_EQ(5, z1[1]);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(Greedy2DTest, find_normal_test)
//{
//    int rval = 0;
//
//    double d1[] = {1.0, 0.0};
//    double d2[] = {1.0, 1.0};
//
//    double x1[] = {1.0,  1.0};
//    double x2[] = {1.0, -1.0};
//
//    double n[2];
//
//    rval = find_normal(d1, x1, n);
//    abort_if(rval, "find_normal failed");
//    EXPECT_DOUBLE_EQ(0, n[0]);
//    EXPECT_DOUBLE_EQ(1, n[1]);
//
//    rval = find_normal(d1, x2, n);
//    abort_if(rval, "find_normal failed");
//    EXPECT_DOUBLE_EQ(0, n[0]);
//    EXPECT_DOUBLE_EQ(-1, n[1]);
//
//    rval = find_normal(d2, x1, n);
//    abort_if(rval, "find_normal failed");
//    EXPECT_DOUBLE_EQ(0, n[0]);
//    EXPECT_DOUBLE_EQ(0, n[1]);
//
//    rval = find_normal(d2, x2, n);
//    abort_if(rval, "find_normal failed");
//    EXPECT_DOUBLE_EQ(1, n[0]);
//    EXPECT_DOUBLE_EQ(-1, n[1]);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(Greedy2DTest, check_rays_parallel)
//{
//    double r1[] = {1.0, 2.0};
//    double r2[] = {2.0, 4.0};
//    double r3[] = {-2.0, -4.0};
//    double r4[] = {0.0, 0.0};
//
//    int match;
//    double scale;
//
//    check_rays_parallel(r1, r2, &match, &scale);
//    EXPECT_TRUE(match);
//    EXPECT_DOUBLE_EQ(scale, 0.5);
//
//    check_rays_parallel(r1, r1, &match, &scale);
//    EXPECT_TRUE(match);
//    EXPECT_DOUBLE_EQ(scale, 1.0);
//
//    check_rays_parallel(r2, r3, &match, &scale);
//    EXPECT_FALSE(match);
//
//    check_rays_parallel(r1, r4, &match, &scale);
//    EXPECT_FALSE(match);
//}
//
//TEST(Greedy2DTest, find_ray)
//{
//    double rays[] = {1.0, 2.0, -1.0, 0.0, -5.0, -5.0};
//    int nrays = 3;
//
//    double r1[] = {1.0, 2.0};
//    double r2[] = {-2.0, 0.0};
//    double r3[] = {-1.0, -1.0};
//    double r4[] = {7.0, 1.0};
//    double r5[] = {-1.0, -2.0};
//
//    int found, index;
//    double scale;
//
//    find_ray(rays, nrays, r1, &found, &scale, &index);
//    EXPECT_TRUE(found);
//    EXPECT_EQ(index, 0);
//    EXPECT_DOUBLE_EQ(scale, 1.0);
//
//    find_ray(rays, nrays, r2, &found, &scale, &index);
//    EXPECT_TRUE(found);
//    EXPECT_EQ(index, 1);
//    EXPECT_DOUBLE_EQ(scale, 2.0);
//
//    find_ray(rays, nrays, r3, &found, &scale, &index);
//    EXPECT_TRUE(found);
//    EXPECT_EQ(index, 2);
//    EXPECT_DOUBLE_EQ(scale, 1 / 5.0);
//
//    find_ray(rays, nrays, r4, &found, &scale, &index);
//    EXPECT_FALSE(found);
//
//    find_ray(rays, nrays, r5, &found, &scale, &index);
//    EXPECT_FALSE(found);
//}
//
//TEST(Greedy2DTest, extract_rays_from_two_sparse_rows_test)
//{
//    int rval = 0;
//
//    char ctypes[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//    double pi1[] = {-3.0, 1.0, -2.0, -1.5, -3.0};
//    int indices1[] = {2, 3, 6, 7, 9};
//
//    double pi2[] = {1.0, -2.0, -3.0, -1.0};
//    int indices2[] = {0, 2, 5, 7};
//
//    struct Row row1 = {.nz = 5, .head = 3, .pi_zero = 0, .pi = pi1, .indices = indices1};
//    struct Row row2 = {.nz = 4, .head = 0, .pi_zero = 0, .pi = pi2, .indices = indices2};
//
//    int nrays;
//    int nz;
//
//    double rays[100];
//    int indices[100];
//    double ray_scale[100];
//    int variable_to_ray[100];
//
//    rval = extract_rays_from_two_sparse_rows(&row1, &row2, ctypes, rays, &nrays,
//            variable_to_ray, ray_scale, indices, &nz);
//    abort_if(rval, "extract_rays_from_two_sparse_rows failed");
//
//    EXPECT_EQ(nrays, 3);
//    EXPECT_EQ(nz, 5);
//
//    EXPECT_DOUBLE_EQ(rays[0], 3.0);
//    EXPECT_DOUBLE_EQ(rays[1], 2.0);
//    EXPECT_DOUBLE_EQ(rays[2], 0.0);
//    EXPECT_DOUBLE_EQ(rays[3], 3.0);
//    EXPECT_DOUBLE_EQ(rays[4], 2.0);
//    EXPECT_DOUBLE_EQ(rays[5], 0.0);
//
//    EXPECT_EQ(indices[0], 2);
//    EXPECT_EQ(indices[1], 5);
//    EXPECT_EQ(indices[2], 6);
//    EXPECT_EQ(indices[3], 7);
//    EXPECT_EQ(indices[4], 9);
//
//    EXPECT_EQ(variable_to_ray[0], 0);
//    EXPECT_EQ(variable_to_ray[1], 1);
//    EXPECT_EQ(variable_to_ray[2], 2);
//    EXPECT_EQ(variable_to_ray[3], 0);
//    EXPECT_EQ(variable_to_ray[4], 2);
//
//    EXPECT_DOUBLE_EQ(ray_scale[0], 1.0);
//    EXPECT_DOUBLE_EQ(ray_scale[1], 1.0);
//    EXPECT_DOUBLE_EQ(ray_scale[2], 1.0);
//    EXPECT_DOUBLE_EQ(ray_scale[3], 0.5);
//    EXPECT_DOUBLE_EQ(ray_scale[4], 1.5);
//
//    CLEANUP:
//    if (rval) FAIL();
//}
//
