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
#include <math.h>
#include <multirow/util.h>
#include <multirow/geometry.h>
}

TEST(GeometryTest, interpolate_coefficients_1)
{
    int rval = 0;

    double r1[] = {1.0, 2.0};
    double r2[] = {3.0, 1.0};
    double p[] = {2.0, 2.0};

    double lambda;

    rval = scale_vector_to_line(r1, r2, p, &lambda);
    abort_if(rval, "scale_vector_to_line failed");

    EXPECT_DOUBLE_EQ(5 / 6.0, lambda);

    CLEANUP:
    if (rval) FAIL();
}

TEST(GeometryTest, interpolate_coefficients_2)
{
    int rval = 0;

    double r1[] = {1.0, 3.0};
    double r2[] = {1e90, 1e90};

    double p1[] = {2.0, 3.0};
    double p2[] = {1.0, 1.0};
    double p3[] = {1.0, 3.0};

    double lambda1;
    double lambda2;
    double lambda3;

    rval = scale_vector_to_line(r1, r2, p1, &lambda1);
    abort_if(rval, "scale_vector_to_line failed");

    rval = scale_vector_to_line(r1, r2, p2, &lambda2);
    abort_if(rval, "scale_vector_to_line failed");

    rval = scale_vector_to_line(r1, r2, p3, &lambda3);
    abort_if(rval, "scale_vector_to_line failed");

    EXPECT_DOUBLE_EQ(2.0, lambda1);
    EXPECT_DOUBLE_EQ(INFINITY, lambda2);
    EXPECT_DOUBLE_EQ(1.0, lambda3);

    rval = scale_vector_to_line(r2, r1, p1, &lambda1);
    abort_if(rval, "scale_vector_to_line failed");

    rval = scale_vector_to_line(r2, r1, p2, &lambda2);
    abort_if(rval, "scale_vector_to_line failed");

    rval = scale_vector_to_line(r2, r1, p3, &lambda3);
    abort_if(rval, "scale_vector_to_line failed");

    EXPECT_DOUBLE_EQ(2.0, lambda1);
    EXPECT_DOUBLE_EQ(INFINITY, lambda2);
    EXPECT_DOUBLE_EQ(1.0, lambda3);

    CLEANUP:
    if (rval) FAIL();
}

TEST(GeometryTest, interpolate_coefficients_3)
{
    int rval = 0;

    double r1[] = {-1.0, 1.0};
    double r2[] = {-1.0, -1.0};

    double p1[] = {-1.0, 0.0};
    double lambda1;

    rval = scale_vector_to_line(r1, r2, p1, &lambda1);
    abort_if(rval, "scale_vector_to_line failed");

    EXPECT_DOUBLE_EQ(1.0, lambda1);

    CLEANUP:
    if (rval) FAIL();
}

TEST(GeometryTest, chull_2d_test)
{
    int rval = 0;

    int npoints = 9;
    double points[] =
    {
        0.0, 0.0,
        0.0, 1.0,
        0.0, 2.0,
        1.0, 0.0,
        1.0, 1.0,
        1.0, 2.0,
        2.0, 0.0,
        2.0, 1.0,
        2.0, 2.0
    };

    double vertices[100];
    int nvertices;

    rval = chull_2d(points, npoints, vertices, &nvertices);
    abort_if(rval, "chull_2d failed");

    for(int i=0; i<nvertices; i++)
        log_debug("  %.2lf %.2lf\n", vertices[2*i], vertices[2*i+1]);

    CLEANUP:
    if(rval) FAIL();
}

//TEST(GeometryTest, sort_rays_test)
//{
//    double rays[] = {1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, -1.0};
//    int var_to_ray[] = {0, 0, 1, 2, 3};
//    double bounds[] = {1.0, 2.0, 3.0, 4.0};
//
//    sort_rays(rays, 4, bounds, var_to_ray, 5);
//
//    EXPECT_EQ(rays[0], -1.0);
//    EXPECT_EQ(rays[1], -1.0);
//
//    EXPECT_EQ(rays[2], -1.0);
//    EXPECT_EQ(rays[3], 1.0);
//
//    EXPECT_EQ(rays[4], 1.0);
//    EXPECT_EQ(rays[5], 1.0);
//
//    EXPECT_EQ(rays[6], 1.0);
//    EXPECT_EQ(rays[7], -1.0);
//
//    EXPECT_EQ(var_to_ray[0], 2);
//    EXPECT_EQ(var_to_ray[1], 2);
//    EXPECT_EQ(var_to_ray[2], 1);
//    EXPECT_EQ(var_to_ray[3], 3);
//    EXPECT_EQ(var_to_ray[4], 0);
//
//    EXPECT_EQ(bounds[0], 4.0);
//    EXPECT_EQ(bounds[1], 2.0);
//    EXPECT_EQ(bounds[2], 1.0);
//    EXPECT_EQ(bounds[3], 3.0);
//}
//
//TEST(GeometryTest, scale_rays)
//{
//    int rval = 0;
//
//    double rays[] = {1.0, 1.0, 2.0, 1.0};
//    int nrays = 2;
//    double scale[] = {5.0, 2.0};
//    
//    rval = scale_rays(rays, nrays, scale);
//    abort_if(rval, "scale_rays failed");
//
//    EXPECT_DOUBLE_EQ(rays[0], 5.0);
//    EXPECT_DOUBLE_EQ(rays[1], 5.0);
//    EXPECT_DOUBLE_EQ(rays[2], 4.0);
//    EXPECT_DOUBLE_EQ(rays[3], 2.0);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, scale_cone_to_point_test)
//{
//    int rval = 0;
//
//    double r1[2] = { 1.0, 0.0 };
//    double r2[2] = { 0.0, 1.0 };
//
//    double p1[2] = { 1.0, 1.0 };
//    double p2[2] = { 1.0, 0.0 };
//    double p3[2] = { -1.0, 1.0 };
//    double p4[2] = { -1.0, 0.0 };
//    
//    double lambda;
//
//    rval = scale_cone_to_point(r1, r2, p1, &lambda);
//    abort_if(rval, "scale_cone_to_point failed");
//    EXPECT_DOUBLE_EQ(2.0, lambda);
//
//    rval = scale_cone_to_point(r1, r2, p2, &lambda);
//    abort_if(rval, "scale_cone_to_point failed");
//    EXPECT_DOUBLE_EQ(1.0, lambda);
//
//    rval = scale_cone_to_point(r1, r2, p3, &lambda);
//    abort_if(rval, "scale_cone_to_point failed");
//    EXPECT_DOUBLE_EQ(0.0, lambda);
//
//    rval = scale_cone_to_point(r1, r2, p4, &lambda);
//    abort_if(rval, "scale_cone_to_point failed");
//    EXPECT_DOUBLE_EQ(-1.0, lambda);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, shear_cone_to_point_test)
//{
//    int rval = 0;
//
//    double r1[2] = { 0.0, 1.0 };
//    double r2[2] = { 1.0, 0.0 };
//
//    double p1[2] = { 1.0, 0.5 };
//    double p2[2] = { 1.5, 0.5 };
//    double p3[2] = { 0.5, 0.0 };
//    double p4[2] = { 1.0, 1.0 };
//    double p5[2] = { 1.0, 2.0 };
//    
//    double lambda;
//
//    rval = shear_cone_to_point(r1, r2, p1, &lambda);
//    abort_if(rval, "shear_cone_to_point failed");
//    EXPECT_DOUBLE_EQ(2.0, lambda);
//
//    rval = shear_cone_to_point(r1, r2, p2, &lambda);
//    abort_if(rval, "shear_cone_to_point failed");
//    EXPECT_DOUBLE_EQ(3.0, lambda);
//
//    rval = shear_cone_to_point(r1, r2, p3, &lambda);
//    abort_if(rval, "shear_cone_to_point failed");
//    EXPECT_DOUBLE_EQ(0.5, lambda);
//
//    rval = shear_cone_to_point(r1, r2, p4, &lambda);
//    abort_if(rval, "shear_cone_to_point failed");
//    EXPECT_DOUBLE_EQ(INFINITY, lambda);
//
//    rval = shear_cone_to_point(r1, r2, p5, &lambda);
//    abort_if(rval, "shear_cone_to_point failed");
//    EXPECT_DOUBLE_EQ(-1.0, lambda);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, compute_bound_for_point_test_1)
//{
//    int rval = 0;
//
//    int nrays = 4;
//    double f[] = { 0.75, 0.75 };
//    double rays[] = { 0.0, 1.0, 1.0, 0.0, 0.0, -1.0, -1.0, 0.0 };
//    double bounds[] = { INFINITY, INFINITY, INFINITY, INFINITY };
//
//    double p1[] = { 0.0, 0.0 };
//    double p2[] = { 0.0, 1.0 };
//    double p3[] = { 1.0, 0.0 };
//    double p4[] = { 1.0, 1.0 };
//
//    double epsilon;
//    double v1[2], v2[2];
//    int i1, i2;
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p1, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(1.5, epsilon);
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p2, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(1.0, epsilon);
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p3, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(1.0, epsilon);
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p4, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(0.5, epsilon);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, compute_bound_for_point_test_2)
//{
//    int rval = 0;
//
//    int nrays = 4;
//
//    double f[] = { 0.75, 0.75 };
//    double rays[] = { 0.0, 1.0, 1.0, 0.0, 0.0, -1.0, -1.0, 0.0 };
//    double bounds[] = { 0.5, 0.5, INFINITY, INFINITY };
//
//    double p1[] = { 0.0, 0.0 };
//    double p2[] = { 0.0, 1.0 };
//    double p3[] = { 1.0, 0.0 };
//    double p4[] = { 1.0, 1.0 };
//
//    double epsilon;
//    double v1[2], v2[2];
//    int i1, i2;
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p1, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(1.5, epsilon);
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p2, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(1.5, epsilon);
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p3, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(1.5, epsilon);
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p4, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(INFINITY, epsilon);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, compute_bound_for_point_test_3)
//{
//    int rval = 0;
//
//    int nrays = 4;
//
//    double f[] = { 0.75, 0.75 };
//    double rays[] = { 0.0, 1.0, 1.0, 0.0, 0.0, -1.0, -1.0, 0.0 };
//    double bounds[] = { 0.5, 0.5, 1.5, 1.5 };
//
//    double p1[] = { 0.0, 0.0 };
//    double p2[] = { 0.0, 1.0 };
//    double p3[] = { 1.0, 0.0 };
//    double p4[] = { 1.0, 1.0 };
//
//    double epsilon;
//    double v1[2], v2[2];
//    int i1, i2;
//
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p1, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(INFINITY, epsilon);
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p2, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(INFINITY, epsilon);
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p3, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(INFINITY, epsilon);
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p4, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(INFINITY, epsilon);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, compute_bound_for_point_test_4)
//{
//    int rval = 0;
//
//    int nrays = 4;
//
//    double f[] = { 0.25, 0.5 };
//    double rays[] = { 1.0, 1.0, 0.0, -1.0, -1.0, 1.0 };
//    double bounds[] = { 0.5, INFINITY, 0.5, INFINITY };
//
//    double p1[] = { 0.0, 0.0 };
//
//    double epsilon;
//    double v1[2], v2[2];
//    int i1, i2;
//
//    rval = compute_bound_for_point(rays, bounds, nrays, f, p1, &epsilon, v1, v2, &i1, &i2);
//    abort_if(rval, "compute_bound_for_point failed");
//    EXPECT_DOUBLE_EQ(1.5, epsilon);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, compute_bounds_2d_test_1)
//{
//    int rval = 0;
//
//    int nrays = 4;
//    double f[] = { 0.75, 0.75 };
//    double rays[] = { 0.0, 1.0, 1.0, 0.0, 0.0, -1.0, -1.0, 0.0 };
//
//    double bounds[4];
//
//    rval = compute_bounds_2d(rays, nrays, f, bounds);
//    abort_if(rval, "compute_bounds_2d failed");
//
//    EXPECT_DOUBLE_EQ(0.5, bounds[0]);
//    EXPECT_DOUBLE_EQ(0.5, bounds[1]);
//    EXPECT_DOUBLE_EQ(1.5, bounds[2]);
//    EXPECT_DOUBLE_EQ(1.5, bounds[3]);
//
//CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, compute_bounds_2d_test_2)
//{
//    int rval = 0;
//
//    int nrays = 3;
//    double f[] = { 0.25, 0.5 };
//    double rays[] = { 1.0, 1.0, 0.0, -1.0, -1.0, 1.0 };
//
//    double bounds[3];
//
//    rval = compute_bounds_2d(rays, nrays, f, bounds);
//    abort_if(rval, "compute_bounds_2d failed");
//
//    EXPECT_DOUBLE_EQ(0.5, bounds[0]);
//    EXPECT_DOUBLE_EQ(1.5, bounds[1]);
//    EXPECT_DOUBLE_EQ(0.5, bounds[2]);
//
//CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, compute_bounds_2d_test_3)
//{
//    int rval = 0;
//
//    int nrays = 5;
//    double f[] = {1 / 4.0, 3 / 4.0};
//    double rays[] = {-2 / 5.0, 5 / 7.0, 0.0, 1.0, 1.0, 1.0, 4 / 5.0, -2 / 3.0,
//            -1.0, 0.0};
//
//    double bounds[5];
//
//    rval = compute_bounds_2d(rays, nrays, f, bounds);
//    abort_if(rval, "compute_bounds_2d failed");
//
//    EXPECT_DOUBLE_EQ(23 / 50.0, bounds[0]);
//    EXPECT_DOUBLE_EQ(23 / 42.0, bounds[1]);
//    EXPECT_DOUBLE_EQ(9 / 11.0, bounds[2]);
//    EXPECT_DOUBLE_EQ(9 / 11.0, bounds[3]);
//    EXPECT_DOUBLE_EQ(23 / 50.0, bounds[4]);
//
//    CLEANUP:
//    if (rval) FAIL();
//}
//
//
//TEST(GeometryTest, compute_bounds_2d_test_4)
//{
//    int rval = 0;
//
//    int nrays = 5;
//    double f[] = {1 / 2.0, 1 / 2.0};
//    double rays[] = {-1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, -1.0, -1.0, -1.0 };
//
//    double bounds[5];
//
//    rval = compute_bounds_2d(rays, nrays, f, bounds);
//    abort_if(rval, "compute_bounds_2d failed");
//
//    EXPECT_DOUBLE_EQ(0.5, bounds[0]);
//    EXPECT_DOUBLE_EQ(0.5, bounds[1]);
//    EXPECT_EQ(GREEDY_BIG_E, bounds[2]);
//    EXPECT_DOUBLE_EQ(0.5, bounds[3]);
//    EXPECT_DOUBLE_EQ(0.5, bounds[4]);
//
//    CLEANUP:
//    if (rval) FAIL();
//}
//
//TEST(GeometryTest, compute_bounds_2d_test_5)
//{
//    int rval = 0;
//
//    int nrays = 6;
//    double f[] = {1 / 2.0, 1 / 2.0};
//    double rays[] = {-1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, -1.0, -1.0, -1.0 };
//
//    double bounds[6];
//
//    rval = compute_bounds_2d(rays, nrays, f, bounds);
//    abort_if(rval, "compute_bounds_2d failed");
//
//    EXPECT_DOUBLE_EQ(0.5, bounds[0]);
//    EXPECT_DOUBLE_EQ(1.0, bounds[1]);
//    EXPECT_DOUBLE_EQ(0.5, bounds[2]);
//    EXPECT_DOUBLE_EQ(1.0, bounds[3]);
//    EXPECT_DOUBLE_EQ(0.5, bounds[4]);
//    EXPECT_DOUBLE_EQ(0.5, bounds[5]);
//
//    CLEANUP:
//    if (rval) FAIL();
//}
//
////TEST(GeometryTest, next_lattice_point_test)
////{
////    int rval = 0;
////
////    LatticeSequence seq;
////    lattice_sequence_init(&seq);
////
////    double lb[] = { -9.0, -1.0 };
////    double ub[] = { 3.0, 4.0 };
////
////
////    while(!seq.eol)
////    {
////        printf("%3d %3d\n", seq.i, seq.j);
////        next_lattice_point(&seq, lb, ub);
////    }
////
////    CLEANUP:
////    if(rval) FAIL();
////}
//
//
//TEST(GeometryTest, get_plane_test)
//{
//    int rval = 0;
//
//    double v0[] = { -1.0, 1.0 };
//    double v1[] = { 1.0, 1.0 };
//    double v2[] = { 0.0, 1.0 };
//    double v3[] = { 1.0, 0.0 };
//
//    double plane[2];
//
//    rval = get_plane(v0, v1, plane);
//    abort_if(rval, "get_plane failed");
//    EXPECT_DOUBLE_EQ(0.0, plane[0]);
//    EXPECT_DOUBLE_EQ(1.0, plane[1]);
//
//    rval = get_plane(v1, v0, plane);
//    abort_if(rval, "get_plane failed");
//    EXPECT_DOUBLE_EQ(0.0, plane[0]);
//    EXPECT_DOUBLE_EQ(1.0, plane[1]);
//
//    rval = get_plane(v2, v3, plane);
//    abort_if(rval, "get_plane failed");
//    EXPECT_DOUBLE_EQ(1.0, plane[0]);
//    EXPECT_DOUBLE_EQ(1.0, plane[1]);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, generate_split_test)
//{
//    int rval = 0;
//
//    double f[2] = { 0.5, 0.0 };
//    double d[2] = { 1.0, 1.0 };
//
//    double pi[2], pi_zero;
//
//    rval = generate_split(f, d, pi, &pi_zero, 1000);
//    abort_if(rval, "generate_split failed");
//
//    EXPECT_DOUBLE_EQ(0.0, pi_zero);
//    EXPECT_DOUBLE_EQ(1.0, pi[0]);
//    EXPECT_DOUBLE_EQ(-1.0, pi[1]);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, generate_split_test_2)
//{
//    int rval = 0;
//
//    double f[2] = { 0.5, 0.25 };
//    double d[2] = { 9/15.0, 13/21.0 };
//
//    double pi[2], pi_zero;
//
//    rval = generate_split(f, d, pi, &pi_zero, 1000);
//    abort_if(rval, "generate_split failed");
//
//    EXPECT_DOUBLE_EQ(16.0, pi_zero);
//    EXPECT_DOUBLE_EQ(65.0, pi[0]);
//    EXPECT_DOUBLE_EQ(-63.0, pi[1]);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, generate_split_test_3)
//{
//    int rval = 0;
//
//    double f[2] = { 0.0, 0.5 };
//    double d[2] = { -1.0, 0.0 };
//
//    double pi[2], pi_zero;
//
//    rval = generate_split(f, d, pi, &pi_zero, 1000);
//    abort_if(rval, "generate_split failed");
//
//    EXPECT_DOUBLE_EQ(0.0, pi_zero);
//    EXPECT_DOUBLE_EQ(0.0, pi[0]);
//    EXPECT_DOUBLE_EQ(1.0, pi[1]);
//
//    CLEANUP:
//    if(rval) FAIL();
//}
//
//TEST(GeometryTest, generate_split_test_4)
//{
//    int rval = 0;
//
//    double f[2] = { 0.5, 0.5 };
//    double d[2] = { 0.0, 1.70 };
//
//    double pi[2], pi_zero;
//
//    rval = generate_split(f, d, pi, &pi_zero, 1000);
//    abort_if(rval, "generate_split failed");
//
//    EXPECT_DOUBLE_EQ(0.0, pi_zero);
//    EXPECT_DOUBLE_EQ(1.0, pi[0]);
//    EXPECT_DOUBLE_EQ(0.0, pi[1]);
//
//    CLEANUP:
//    if(rval) FAIL();
//}


//TEST(GeometryTest, find_linear_combination)
//{
//    double p[] = {1.0, 1.0};
//    double r1[] = {1 / 3.0, 1 / 5.0};
//    double r2[] = {-1.0, 2 / 3.0};
//    double lambda1, lambda2;
//
//    find_linear_combination(r1, r2, p, &lambda1, &lambda2);
//
//    EXPECT_NEAR(3.94736842105263, lambda1, 1e-12);
//    EXPECT_NEAR(0.315789473684211, lambda2, 1e-12);
//}
