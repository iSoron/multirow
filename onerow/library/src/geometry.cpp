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

#include <onerow/geometry.hpp>

Point::Point() :
		x(0), y(0)
{

}

Point::Point(rational _x, rational _y) :
		x(_x), y(_y)
{

}

Point Point::operator+(const Point& p) const
{
	return Point(x + p.x, y + p.y);
}

Point Point::operator-(const Point& p) const
{
	return Point(x - p.x, y - p.y);
}

rational Point::operator*(const Point& p) const
{
	return x * p.x + y * p.y;
}

Point Point::operator*(rational scale) const
{
	return Point(scale * x, scale * y);
}

Line::Line()
{

}

Line::Line(Point _p1, Point _p2) :
		p1(_p1), p2(_p2)
{

}

Line::Line(rational x1, rational y1, rational x2, rational y2) :
		p1(Point(x1, y1)), p2(Point(x2, y2))
{

}

std::ostream& operator<<(std::ostream& os, const Point &p)
{
	os << "(" << p.x << "," << p.y << ")";
	return os;
}

std::ostream& operator<<(std::ostream& os, const Line &l)
{
	os << l.p1 << "--" << l.p2;
	return os;
}
