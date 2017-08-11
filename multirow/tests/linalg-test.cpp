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

extern "C" {
#include <multirow/util.h>
#include <multirow/linalg.h>
}

const double E = 1e-6;

TEST(LinAlgTest, scale_test)
{
    double x[] = {1.0, 2.0, 3.0};
    LINALG_scale(3, x, 2.0);
    EXPECT_NEAR(x[0], 2.0, E);
    EXPECT_NEAR(x[1], 4.0, E);
    EXPECT_NEAR(x[2], 6.0, E);
}

TEST(LinAlgTest, dot_test)
{
    double x[] = { 1.0, 2.0, 3.0 };
    double y[] = { 3.0, 4.0, 5.0 };
    double dot = LINALG_dot(3, x, y);
    EXPECT_NEAR(dot, 26.0, E);
}

TEST(LinAlgTest, norm_test)
{
    double x[] = { 1.0,  2.0, -3.0 };
    double y[] = { 3.0, -4.0, -5.0 };
    double x_norm = LINALG_norm(3, x);
    double y_norm = LINALG_norm(3, y);
    EXPECT_NEAR(x_norm, 6.0, E);
    EXPECT_NEAR(y_norm, 12.0, E);
}

TEST(LinAlgTest, solve_test)
{
    int rval = 0;

    double A[] = {
            2.0, 1.0,  3.0,
            2.0, 6.0,  8.0,
            6.0, 8.0, 18.0,
    };
    double b[] = { 1.0, 3.0, 5.0 };
    double x[] = { 0.0, 0.0, 0.0 };

    rval = LINALG_solve(3, A, b, x);
    abort_if(rval, "LINALG_solve failed");

    // Should compute x correctly
    EXPECT_NEAR(x[0], 0.3, E);
    EXPECT_NEAR(x[1], 0.4, E);
    EXPECT_NEAR(x[2], 0.0, E);

    // Should not modify A and b
    EXPECT_EQ(A[0], 2.0);
    EXPECT_EQ(A[1], 1.0);
    EXPECT_EQ(A[2], 3.0);
    EXPECT_EQ(b[0], 1.0);
    EXPECT_EQ(b[1], 3.0);
    EXPECT_EQ(b[2], 5.0);

CLEANUP:
    if(rval) FAIL();
}