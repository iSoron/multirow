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
#include <onerow/knapsack2.hpp>
typedef q::mpq rational;

using std::cout;
using std::endl;

TEST(Knapsack2Test, test_1)
{
	Knapsack2 k(rational(5,7), rational(3,5));

//	cout << "list:\n";
//	for(auto v : k.list)
//		cout << v.lower << " "  << v.upper << " " << v.side << " " << v.opposed << endl;

//	cout << "has ray? " << k.has_ray << endl;
}

TEST(Knapsack2Test, test_2)
{
	Knapsack2 k(rational(1,2), rational(1,2));

//	cout << "left:\n";
//	each(k.left, v)
//		cout << v->active << v->opposed << endl;
//
//	cout << "right:\n";
//	each(k.right, v)
//		cout << v->active << v->opposed << endl;
//
//	cout << "has ray? " << k.has_ray << endl;
}
