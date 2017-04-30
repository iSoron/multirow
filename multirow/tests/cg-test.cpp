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

#define TEST_SOURCE

#include <gtest/gtest.h>
#include <math.h>

extern "C" {
#include <multirow/cg.h>
#include "../src/cg.c"
}

int BOOST_VAR = -1;
double BOOST_FACTOR = 1.0;

TEST(CGTest, next_combination_test_1)
{
    int rval = 0;

    int v[] = { 1, 0 };
    int n = 10;
    int k = 2;

    int finished;
    int count = 0;

    do
    {
        count++;
        next_combination(n, k, v, 0, &finished);
    }
    while(!finished);

    EXPECT_EQ(count, 45);

    CLEANUP:
    if(rval) FAIL();
}

TEST(CGTest, next_combination_test_2)
{
    int rval = 0;

    int v[] = { 2, 1, 0 };
    int n = 5;
    int k = 3;

    int finished;
    int count = 0;

    do
    {
        count++;
        next_combination(n, k, v, 0, &finished);
    }
    while(!finished);

    EXPECT_EQ(count, 10);

CLEANUP:
    if(rval) FAIL();
}

TEST(CGTest, check_rays_parallel_test)
{
    int rval = 0;

    double r1[] = { 1.0, 1.0, 1.0 };
    double r2[] = { 2.0, 2.0, 2.0 };
    double r3[] = { 0.0, 0.0, 0.0 };
    double r4[] = { -1.0, -1.0, -1.0 };

    int match;
    double scale;

    rval = check_rays_parallel(3, r1, r2, &match, &scale);
    abort_if(rval, "check_rays_parallel failed");
    EXPECT_TRUE(match);
    EXPECT_DOUBLE_EQ(scale, 0.5);

    rval = check_rays_parallel(3, r1, r3, &match, &scale);
    abort_if(rval, "check_rays_parallel failed");
    EXPECT_FALSE(match);

    rval = check_rays_parallel(3, r1, r4, &match, &scale);
    abort_if(rval, "check_rays_parallel failed");
    EXPECT_FALSE(match);

CLEANUP:
    if(rval) FAIL();

}

TEST(CGTest, extract_rays_from_rows_test)
{
    int rval = 0;

    char column_types[16] = {0};
    double pi1[] = { 1.0, 1.0, 1.0, 2.0, 1.0 };
    int indices1[] = { 1, 7, 8, 12, 14 };
    struct Row row1 =
    {
        .nz = 5,
        .head = 1,
        .pi_zero = 1.0,
        .pi = pi1,
        .indices = indices1
    };
    
    double pi2[] = { 1.0, 1.0, 2.0, 1.0, 1.0, 1.0 };
    int indices2[] = { 5, 8, 12, 13, 14, 16 };
    struct Row row2 =
    {
        .nz = 6,
        .head = 13,
        .pi_zero = 1.0,
        .pi = pi2,
        .indices = indices2
    };

    double pi3[] = { 1.0, 1.0, 1.0, 1.0 };
    int indices3[] = { 3, 7, 10, 16 };
    struct Row row3 =
    {
        .nz = 4,
        .head = 7,
        .pi_zero = 1.0,
        .pi = pi3,
        .indices = indices3
    };

    struct Row *rows[] = { &row1, &row2, &row3 };
    struct Tableau tableau = {3, rows, column_types};

    int indices[1000];
    int variable_to_ray[1000];
    double ray_scale[1000];

    struct RayList rays;
    LFREE_init_ray_list(&rays, 3, 1000);
    struct RayMap map = {rays, variable_to_ray, ray_scale, indices, 0};

    rval = CG_extract_rays_from_tableau(&tableau, &map);
    abort_if(rval, "CG_extract_rays_from_rows failed");

    EXPECT_EQ(rays.nrays, 4);
    EXPECT_EQ(map.nvars, 7);

    EXPECT_DOUBLE_EQ(rays.values[0],  0.0);
    EXPECT_DOUBLE_EQ(rays.values[1],  0.0);
    EXPECT_DOUBLE_EQ(rays.values[2], -1.0);

    EXPECT_DOUBLE_EQ(rays.values[3],  0.0);
    EXPECT_DOUBLE_EQ(rays.values[4], -1.0);
    EXPECT_DOUBLE_EQ(rays.values[5],  0.0);

    EXPECT_DOUBLE_EQ(rays.values[6], -2.0);
    EXPECT_DOUBLE_EQ(rays.values[7], -2.0);
    EXPECT_DOUBLE_EQ(rays.values[8],  0.0);

    EXPECT_DOUBLE_EQ(rays.values[ 9],  0.0);
    EXPECT_DOUBLE_EQ(rays.values[10], -1.0);
    EXPECT_DOUBLE_EQ(rays.values[11], -1.0);

    EXPECT_EQ(indices[0], 3);
    EXPECT_EQ(indices[1], 5);
    EXPECT_EQ(indices[2], 8);
    EXPECT_EQ(indices[3], 10);
    EXPECT_EQ(indices[4], 12);
    EXPECT_EQ(indices[5], 14);
    EXPECT_EQ(indices[6], 16);

    EXPECT_EQ(ray_scale[0], 1.0);
    EXPECT_EQ(ray_scale[1], 1.0);
    EXPECT_EQ(ray_scale[2], 0.5);
    EXPECT_EQ(ray_scale[3], 1.0);
    EXPECT_EQ(ray_scale[4], 1.0);
    EXPECT_EQ(ray_scale[5], 0.5);
    EXPECT_EQ(ray_scale[6], 1.0);

    EXPECT_EQ(variable_to_ray[0], 0);
    EXPECT_EQ(variable_to_ray[1], 1);
    EXPECT_EQ(variable_to_ray[2], 2);
    EXPECT_EQ(variable_to_ray[3], 0);
    EXPECT_EQ(variable_to_ray[4], 2);
    EXPECT_EQ(variable_to_ray[5], 2);
    EXPECT_EQ(variable_to_ray[6], 3);

    CLEANUP:
    if(rval) FAIL();
}

TEST(CGTest, boost_variable_test)
{
    int rval = 0;

    int nrows = 3;
    double rays[] = {
         0,  0, -1,
         0, -1,  0,
        -2, -2,  0,
         0, -1, -1
    };

    int nz = 7;
    int indices[] = { 3, 5, 8, 10, 12, 14, 16 };
    int variable_to_ray[] = { 0, 1, 2, 0, 2, 2, 3 };
    double ray_scale[] = { 1, 1, 0.5, 1, 1, 0.5, 1 };

    rval = CG_boost_variable(8, 200, nrows, rays, variable_to_ray, ray_scale, indices, nz);
    abort_if(rval, "boost_variable failed");

    EXPECT_DOUBLE_EQ(rays[0],  0.0);
    EXPECT_DOUBLE_EQ(rays[1],  0.0);
    EXPECT_DOUBLE_EQ(rays[2], -1.0);

    EXPECT_DOUBLE_EQ(rays[3],  0.0);
    EXPECT_DOUBLE_EQ(rays[4], -1.0);
    EXPECT_DOUBLE_EQ(rays[5],  0.0);

    EXPECT_DOUBLE_EQ(rays[6], -400.0);
    EXPECT_DOUBLE_EQ(rays[7], -400.0);
    EXPECT_DOUBLE_EQ(rays[8],    0.0);

    EXPECT_DOUBLE_EQ(rays[ 9],  0.0);
    EXPECT_DOUBLE_EQ(rays[10], -1.0);
    EXPECT_DOUBLE_EQ(rays[11], -1.0);

    EXPECT_EQ(indices[0], 3);
    EXPECT_EQ(indices[1], 5);
    EXPECT_EQ(indices[2], 8);
    EXPECT_EQ(indices[3], 10);
    EXPECT_EQ(indices[4], 12);
    EXPECT_EQ(indices[5], 14);
    EXPECT_EQ(indices[6], 16);

    EXPECT_EQ(ray_scale[0], 1.0);
    EXPECT_EQ(ray_scale[1], 1.0);
    EXPECT_EQ(ray_scale[2], 0.0025);
    EXPECT_EQ(ray_scale[3], 1.0);
    EXPECT_EQ(ray_scale[4], 0.005);
    EXPECT_EQ(ray_scale[5], 0.0025);
    EXPECT_EQ(ray_scale[6], 1.0);

    EXPECT_EQ(variable_to_ray[0], 0);
    EXPECT_EQ(variable_to_ray[1], 1);
    EXPECT_EQ(variable_to_ray[2], 2);
    EXPECT_EQ(variable_to_ray[3], 0);
    EXPECT_EQ(variable_to_ray[4], 2);
    EXPECT_EQ(variable_to_ray[5], 2);
    EXPECT_EQ(variable_to_ray[6], 3);

CLEANUP:
    if(rval) FAIL();
}
