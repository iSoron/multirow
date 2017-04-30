/* Copyright (c) 2015-2017 Alinson Xavier
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
#include "../src/infinity.c"
}

TEST(InfinityTest, cmp_ray_angle_test)
{
    double r0[] = { 1.0, 0.0 };
    double r1[] = { 2.0, 0.0 };
    double r2[] = { 1.0, 1.0 };
    double r3[] = { -1.0, 0.0 };
    double r4[] = { 1.0, -1.0 };

    SortPair sp0 = { 0, &r0 };
    SortPair sp1 = { 1, &r1 };
    SortPair sp2 = { 2, &r2 };
    SortPair sp3 = { 3, &r3 };
    SortPair sp4 = { 4, &r4 };

    EXPECT_EQ(_qsort_cmp_rays_angle(&sp0, &sp1),  0);
    EXPECT_EQ(_qsort_cmp_rays_angle(&sp0, &sp2),  1);
    EXPECT_EQ(_qsort_cmp_rays_angle(&sp0, &sp3),  1);
    EXPECT_EQ(_qsort_cmp_rays_angle(&sp0, &sp4), -1);

    EXPECT_EQ(_qsort_cmp_rays_angle(&sp2, &sp0), -1);
    EXPECT_EQ(_qsort_cmp_rays_angle(&sp2, &sp1), -1);
    EXPECT_EQ(_qsort_cmp_rays_angle(&sp2, &sp3),  1);
    EXPECT_EQ(_qsort_cmp_rays_angle(&sp2, &sp4), -1);

    EXPECT_EQ(_qsort_cmp_rays_angle(&sp3, &sp0), -1);
}

TEST(InfinityTest, sort_rays_angle_test)
{
    int rval = 0;

    double r0[] = { 1.0, 1.0 };
    double r1[] = { 1.0, 0.0 };
    double r2[] = { 1.0, -1.0 };
    double r3[] = { -1.0, 0.0 };
    double r4[] = { 2.0, 0.0 };

    RayList rays;
    LFREE_init_ray_list(&rays, 2, 5);
    LFREE_push_ray(&rays, r0);
    LFREE_push_ray(&rays, r1);
    LFREE_push_ray(&rays, r2);
    LFREE_push_ray(&rays, r3);
    LFREE_push_ray(&rays, r4);

    rval = sort_rays_by_angle(&rays);
    abort_if(rval, "sort_rays_by_angle failed");

    SortPair sp0, sp1, sp2, sp3, sp4;
    sp0 = { 0, LFREE_get_ray(&rays, 0) };
    sp1 = { 1, LFREE_get_ray(&rays, 1) };
    sp2 = { 2, LFREE_get_ray(&rays, 2) };
    sp3 = { 3, LFREE_get_ray(&rays, 3) };
    sp4 = { 4, LFREE_get_ray(&rays, 4) };

    EXPECT_LE(_qsort_cmp_rays_angle(&sp0, &sp1),  0);
    EXPECT_LE(_qsort_cmp_rays_angle(&sp1, &sp2),  0);
    EXPECT_LE(_qsort_cmp_rays_angle(&sp2, &sp3),  0);
    EXPECT_LE(_qsort_cmp_rays_angle(&sp3, &sp4),  0);

CLEANUP:
    LFREE_free_ray_list(&rays);
    if (rval) FAIL();
}
