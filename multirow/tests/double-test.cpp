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
#include <math.h>

extern "C" {
    #include <multirow/double.h>
    #include <multirow/util.h>
}

TEST(DoubleTest, double_cmp_test)
{
    EXPECT_EQ(0, DOUBLE_cmp(1/3.0, 2/6.0));


    EXPECT_EQ(1, DOUBLE_cmp(5.0, 3.0));
    EXPECT_EQ(0, DOUBLE_cmp(5.0, 5.0));
    EXPECT_EQ(-1, DOUBLE_cmp(5.0, 10.0));

    EXPECT_EQ(1, DOUBLE_cmp(5.1, 5.0));

    EXPECT_EQ(0, DOUBLE_cmp(1000000005, 1000000003));
    EXPECT_EQ(0, DOUBLE_cmp(1000000005, 1000000005));
    EXPECT_EQ(0, DOUBLE_cmp(1000000005, 1000000007));


    EXPECT_EQ(1, DOUBLE_cmp(0.00005, 0.00003));
    EXPECT_EQ(0, DOUBLE_cmp(0.00005, 0.00005));
    EXPECT_EQ(-1, DOUBLE_cmp(0.00005, 0.00007));

    EXPECT_EQ(0, DOUBLE_cmp(0.5 * EPSILON, 0.3 * EPSILON));
    EXPECT_EQ(0, DOUBLE_cmp(0.5 * EPSILON, 0.5 * EPSILON));
    EXPECT_EQ(0, DOUBLE_cmp(0.5 * EPSILON, 0.7 * EPSILON));

    EXPECT_EQ(0, DOUBLE_cmp(0.0, 0.0));
    EXPECT_EQ(0, DOUBLE_cmp(0.0, -0.0));
    EXPECT_EQ(0, DOUBLE_cmp(-0.0, -0.0));
    EXPECT_EQ(0, DOUBLE_cmp(0.1 * EPSILON, 0.0));
    EXPECT_EQ(0, DOUBLE_cmp(0.0, 0.1 * EPSILON));
    EXPECT_EQ(0, DOUBLE_cmp(-0.1 * EPSILON, 0.0));
    EXPECT_EQ(0, DOUBLE_cmp(0.0, -0.1 * EPSILON));

    EXPECT_EQ(0, DOUBLE_cmp(DBL_MAX, DBL_MAX));
    EXPECT_EQ(1, DOUBLE_cmp(DBL_MAX, -DBL_MAX));
    EXPECT_EQ(-1, DOUBLE_cmp(-DBL_MAX, DBL_MAX));
    EXPECT_EQ(1, DOUBLE_cmp(DBL_MAX, DBL_MAX / 2));
    EXPECT_EQ(1, DOUBLE_cmp(DBL_MAX, -DBL_MAX / 2));
    EXPECT_EQ(-1, DOUBLE_cmp(-DBL_MAX, DBL_MAX / 2));

    EXPECT_EQ(0, DOUBLE_cmp(INFINITY, INFINITY));
    EXPECT_EQ(1, DOUBLE_cmp(INFINITY, -INFINITY));
    EXPECT_EQ(-1, DOUBLE_cmp(-INFINITY, INFINITY));
    EXPECT_EQ(0, DOUBLE_cmp(-INFINITY, -INFINITY));
    EXPECT_EQ(1, DOUBLE_cmp(INFINITY, DBL_MAX));
    EXPECT_EQ(-1, DOUBLE_cmp(-INFINITY, DBL_MIN));

    EXPECT_EQ(1, DOUBLE_cmp(NAN, NAN));
    EXPECT_EQ(1, DOUBLE_cmp(NAN, 0.0f));
    EXPECT_EQ(1, DOUBLE_cmp(-0.0f, NAN));
    EXPECT_EQ(1, DOUBLE_cmp(NAN, -0.0f));
    EXPECT_EQ(1, DOUBLE_cmp(0.0f, NAN));
    EXPECT_EQ(1, DOUBLE_cmp(NAN, INFINITY));
    EXPECT_EQ(1, DOUBLE_cmp(INFINITY, NAN));
    EXPECT_EQ(1, DOUBLE_cmp(NAN, -INFINITY));
    EXPECT_EQ(1, DOUBLE_cmp(-INFINITY, NAN));
    EXPECT_EQ(1, DOUBLE_cmp(NAN, DBL_MAX));
    EXPECT_EQ(1, DOUBLE_cmp(DBL_MAX, NAN));
    EXPECT_EQ(1, DOUBLE_cmp(NAN, -DBL_MAX));
    EXPECT_EQ(1, DOUBLE_cmp(-DBL_MAX, NAN));
    EXPECT_EQ(1, DOUBLE_cmp(NAN, DBL_MIN));
    EXPECT_EQ(1, DOUBLE_cmp(DBL_MIN, NAN));
    EXPECT_EQ(1, DOUBLE_cmp(NAN, -DBL_MIN));
    EXPECT_EQ(1, DOUBLE_cmp(-DBL_MIN, NAN));
}

TEST(DoubleTest, double_iszero_test)
{
    EXPECT_TRUE(DOUBLE_iszero(0.0));
    EXPECT_TRUE(DOUBLE_iszero(-0.0));
    EXPECT_TRUE(DOUBLE_iszero(0.1 * EPSILON));
    EXPECT_TRUE(DOUBLE_iszero(-0.1 * EPSILON));
}

TEST(DoubleTest, double_eq_test)
{
    EXPECT_TRUE(DOUBLE_eq(1.0, 1.0));
    EXPECT_FALSE(DOUBLE_eq(1.0, 3.0));
    EXPECT_TRUE(DOUBLE_eq(1/3.0, 2/6.0));
}

TEST(DoubleTest, double_max_test)
{
    EXPECT_EQ(0.1, DOUBLE_max(0.0, 0.1));
    EXPECT_EQ(0.5, DOUBLE_max(0.5, 0.1));
}

TEST(DoubleTest, double_geq_test)
{
    EXPECT_TRUE(DOUBLE_geq(1.0, 0.5));
    EXPECT_TRUE(DOUBLE_geq(0.46000000, 0.45953916));
}

TEST(DoubleTest, double_to_rational_test)
{
    int rval = 0;

    Rational r;
    
    rval = DOUBLE_to_rational(M_PI, 1, r);
    abort_if(rval, "DOUBLE_to failed");
    EXPECT_EQ(3, r->num);
    EXPECT_EQ(1, r->den);
    
    rval = DOUBLE_to_rational(M_PI, 10, r);
    abort_if(rval, "DOUBLE_to failed");
    EXPECT_EQ(22, r->num);
    EXPECT_EQ(7, r->den);
    
    rval = DOUBLE_to_rational(M_PI, 100, r);
    abort_if(rval, "DOUBLE_to failed");
    EXPECT_EQ(311, r->num);
    EXPECT_EQ(99, r->den);

    rval = DOUBLE_to_rational(M_PI, 1000, r);
    abort_if(rval, "DOUBLE_to failed");
    EXPECT_EQ(355, r->num);
    EXPECT_EQ(113, r->den);

    rval = DOUBLE_to_rational(M_PI, 10000, r);
    abort_if(rval, "DOUBLE_to failed");
    EXPECT_EQ(355, r->num);
    EXPECT_EQ(113, r->den);

    rval = DOUBLE_to_rational(-M_PI, 10000, r);
    abort_if(rval, "DOUBLE_to failed");
    EXPECT_EQ(-355, r->num);
    EXPECT_EQ(113, r->den);

CLEANUP:
    if(rval) FAIL();
}

TEST(DoubleTest, gcd_test)
{
    EXPECT_EQ(1, gcd(11, 7));
    EXPECT_EQ(2, gcd(6, 4));
    EXPECT_EQ(2, gcd(4, 6));
    EXPECT_EQ(5, gcd(5, 10));
}

