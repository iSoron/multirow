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

#ifndef GEOMETRY_HPP_
#define GEOMETRY_HPP_

#include<qxx/rational.hpp>
#include<qxx/svec.hpp>
typedef q::mpq rational;
typedef q::svec svec;

/**
 * Models a two-dimensional rational point (x,y).
 */
struct Point {

	rational x;
	rational y;

	Point();
	Point(rational x, rational y);

	Point operator+(const Point& p) const;
	Point operator-(const Point& p) const;
	Point operator*(rational scale) const;
	rational operator*(const Point& p) const;
};

/**
 * Models a two-dimensional rational line, defined by a pair of points p1 and p2.
 */
struct Line {
	Point p1;
	Point p2;

	Line();
	Line(Point p1, Point p2);
	Line(rational x1, rational y1, rational x2, rational y2);
};

std::ostream& operator<<(std::ostream& os, const Line &l);
std::ostream& operator<<(std::ostream& os, const Point &p);

#endif /* GEOMETRY_HPP_ */
