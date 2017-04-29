/*
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

#include <stdexcept>
#include <gtest/gtest.h>
#include <qxx/rational.hpp>
typedef q::mpq rational;

TEST(RationalTest, reduce_test)
{
	EXPECT_EQ(rational(5), rational(5).reduce(100));
	EXPECT_EQ(rational(0), rational(0).reduce(100));
	EXPECT_EQ(rational(1,3), rational(1/3.0).reduce(100));
	EXPECT_EQ(rational(1,2), rational(1,2).reduce(100));
	EXPECT_EQ(rational(-1,2), rational(-1,2).reduce(100));
}

TEST(RationalTest, minus_test)
{
	EXPECT_EQ(rational(-1,3), -rational(1,3));
}

TEST(RationalTest, plus_minus_times_test)
{
	rational a(47,17);
	rational b(113,98);

	EXPECT_EQ(rational(6527,1666), a+b);
	EXPECT_EQ(rational(2685,1666), a-b);
	EXPECT_EQ(rational(-2685,1666), -a+b);
	EXPECT_EQ(rational(5311,1666), a*b);
}

TEST(RationalTest, floor_ceil_frac_test)
{
	EXPECT_EQ(rational(45,17).floor(), rational(2));
	EXPECT_EQ(rational(45,17).ceil(), rational(3));
	EXPECT_EQ(rational(45,17).frac(), rational(11,17));

	EXPECT_EQ(rational(-45,17).floor(), rational(-3));
	EXPECT_EQ(rational(-45,17).ceil(), rational(-2));
	EXPECT_EQ(rational(-45,17).frac(), rational(6,17));

	EXPECT_EQ(rational(1,2).frac(), rational(1,2));
	EXPECT_EQ(rational(-1,2).frac(), rational(1,2));
}

TEST(RationalTest, double_to_rational_test)
{
	EXPECT_FLOAT_EQ(rational(1/3.0).get_double(), 1/3.0);
	EXPECT_FLOAT_EQ(rational(127/15.0).get_double(), 127/15.0);
	EXPECT_FLOAT_EQ(rational(127/183.0).get_double(), 127/183.0);
	EXPECT_FLOAT_EQ(rational(513/3577.0).get_double(), 513/3577.0);
}
