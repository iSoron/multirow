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
#include <cassert>
#include <cmath>
#include <algorithm>
#include <set>
#include <gmp.h>
#include <onerow/wedge_cut_generator.hpp>
#include <onerow/stats.hpp>
#include <onerow/params.hpp>

using std::cout;
using std::endl;
using std::min;

static bool debug = false;

IntersectionCut::IntersectionCut(Point _f, int _n, const Line *l) :
		f(_f), n_faces(_n), pre_lifting_ready(false)
{
	d = new Point[n_faces];
	if (l != 0)
		for (int i = 0; i < n_faces; i++)
			set_face(i, l[i]);
}

IntersectionCut::~IntersectionCut()
{
	delete[] d;
}

void IntersectionCut::set_face(int index, Line line)
{
	if (index < 0 || index >= n_faces)
		throw std::out_of_range("");

	rational x1 = line.p1.x;
	rational y1 = line.p1.y;
	rational x2 = line.p2.x;
	rational y2 = line.p2.y;

	rational h[3] =	{ -y2 + y1, x2 - x1, x1 * y2 - x2 * y1 };
	rational rhs = -(h[2] + h[0] * f.x + h[1] * f.y);

	if (rhs == rational(0))
		throw std::invalid_argument("");

	d[index].x = h[0] / rhs;
	d[index].y = h[1] / rhs;
}

rational IntersectionCut::get_continuous_coefficient(rational rx, rational ry)
{
	#ifdef ENABLE_EXTENDED_STATISTICS
		Stats::add_coefficient(false);
	#endif

	rational max_coeff(-100000, 1);
	for (int i = 0; i < n_faces; i++)
	{
		rational coeff = rx * d[i].x + ry * d[i].y;
		if (coeff > max_coeff)
			max_coeff = coeff;
	}
	return max_coeff;
}

void IntersectionCut::pre_lifting()
{
	rational r0y;
	rational fx = f.x, fy = f.y;
	rational d0x = d[0].x, d0y = d[0].y;
	rational d1x = d[1].x, d1y = d[1].y;
	rational apex_x, apex_y;
	rational w1 = rational(1) + (d0x * fx + d0y * fy);
	rational w2 = rational(1) + (d1x * fx + d1y * fy);

	rational m = (d0x*d1y-d0y*d1x);

	if (m == 0)
	{
		r0x = -d0y/d0x;
		r0y = 1;
	}
	else
	{

		apex_x = (d1y * w1 - d0y * w2) / m;
		apex_y = (-d1x * w1 + d0x * w2) / m;

		r0x = apex_x - fx;
		r0y = apex_y - fy;

		r0x /= r0y;
		r0y = 1;
	}

	assert(d0x < 0);
	assert(d1x > 0);

	p = r0x*d0x + r0y*d0y;
	assert(p == r0x*d1x + r0y*d1y);
	assert(p >= 0);

	this->d_p = p.get_double();
	this->d_d1x = d1x.get_double();
	this->d_d1y = d1y.get_double();
	this->d_d0x = d0x.get_double();
	this->d_d0y = d0y.get_double();
	this->d_r0x = r0x.get_double();

	pre_lifting_ready = true;
}

double IntersectionCut::get_trivial_lifting_coefficient_double(double rx, double ry)
{
	if(!pre_lifting_ready) pre_lifting();

	#ifdef ENABLE_EXTENDED_STATISTICS
		Stats::add_coefficient(true);
	#endif

	double a = 100000;

	unsigned long k2 = 0;
	unsigned long M = 10000;

	while (k2 < M)
	{
		double b = (k2*d_r0x-rx);

		double a1 = d_d1x * (rx + ceil(b)) + d_d1y * k2;
		double a2 = d_d0x * (rx + floor(b)) + d_d0y * k2;

		assert(rx + ceil(b) + 0.0001 >= d_r0x*k2);
		assert(rx + floor(b) - 0.0001 <= d_r0x*k2);

		assert(a1 >= 0);
		assert(a2 >= 0);

		if(a1 < a) a = a1;
		if(a2 < a) a = a2;

		if(fabs(d_p) < ZERO_CUTOFF) break;

		M = ceil(a/d_p);

		k2++;
	}


	#ifdef ENABLE_EXTENDED_STATISTICS
		Stats::add_trivial_lifting_m(k2);
	#endif

	assert(a >= 0);
	return a;
}


rational IntersectionCut::get_trivial_lifting_coefficient(rational rx, rational ry)
{
	if(!pre_lifting_ready) pre_lifting();

	#ifdef ENABLE_EXTENDED_STATISTICS
		Stats::add_coefficient(true);
	#endif

	rational a = rational(100000, 1);

	unsigned long k2 = 0;
	unsigned long M = 10000;
	rational b = -rx;

	while (k2 < M)
	{
		rational r_k2((signed) k2);
		rational a1 = d[1].x*(rx + b.ceil())  + d[1].y * r_k2;
		rational a2 = d[0].x*(rx + b.floor()) + d[0].y * r_k2;

		assert(rx + b.ceil()  >= r0x*r_k2);
		assert(rx + b.floor() <= r0x*r_k2);

		assert(a1 >= 0);
		assert(a2 >= 0);

		if(a1 < a) a = a1;
		if(a2 < a) a = a2;

		if(p == 0) break;
		M = (a/p).ceil().get_double();

		k2++;
		b += r0x;
	}

	return a;
}


WedgeCut::WedgeCut(Point _f, Point left, Point apex, Point right) :
		IntersectionCut(_f, 2)
{
	set_face(0, Line(left, apex));
	set_face(1, Line(apex, right));
}

SplitCut::SplitCut(Point _f, Point left, Point right, Point direction) :
		IntersectionCut(_f, 2)
{
	set_face(0, Line(left, left + direction));
	set_face(1, Line(right, right + direction));
}

WedgeCutGenerator::WedgeCutGenerator(Row &r) :
		SingleRowCutGenerator(r), finished(false),
		f(2), r1(2),
		r1_offset(-1),
		cur_facet(-1),
		n_knapsacks(0)
{
	f[0] = r.c.pi_zero.frac();
	f[1] = 0;
	cut = new Constraint;
	
	eval_next();
}

WedgeCutGenerator::~WedgeCutGenerator()
{
	delete cut;
}

bool WedgeCutGenerator::has_next()
{
	return !finished;
}

Constraint* WedgeCutGenerator::next()
{
	if (!has_next())
		throw std::out_of_range("");
	
	Constraint *old_cut = cut;

	cut = new Constraint;
	eval_next();

	return(old_cut);
}

q::dvec WedgeCutGenerator::intersection(
	q::dvec a, q::dvec b, q::dvec c, q::dvec d)
{
	q::dmat g(2,2);
	
	g.set_col(0, b - a);
	g.set_col(1, c - d);
	
	q::dvec mult = g.inv() * (c - a);
	
	q::dvec x = a + (b - a) * mult[0];
	
	assert(x == c + (d - c) * mult[1]);
	
	return(x);
}

void WedgeCutGenerator::eval_next()
{
	while (true)
	{
		if (0 <= cur_facet && cur_facet < (int) knapsack.list.size() && cur_facet <= MAX_CUT_DEPTH)
			break;

		r1_offset++;

		if (r1_offset >= row.c.pi.nz())
		{
			finished = true;
			return;
		}

		int r1_index = row.c.pi.index(r1_offset);
		if(row.reduced_costs[r1_index] < row.cost_cutoff)
			continue;

		r1[0] = -row.c.pi.value(r1_offset).reduce(REDUCE_FACTOR_R1);
		r1[1] = 1;

		if (r1_index == row.basic_var_index)
			continue;

		if (!row.is_integer[r1_index])
			continue;

		if (r1[0] == 0)
			continue;

		knapsack.clear();
		knapsack.eval(f[0], r1[0]);

		n_knapsacks++;
		cur_facet = 0;
	}

	Line lines[2];
	int side = knapsack.list[cur_facet].side;
	q::dvec &a1 = knapsack.list[cur_facet].lower;
	q::dvec &a2 = knapsack.list[cur_facet].upper;
	q::dvec &o = knapsack.list[cur_facet].opposed;

	if (side == KNAPSACK2_RAY)
	{
		a2 = o;

		q::dvec b1 = a1 + r1;
		q::dvec b2 = a2 + r1;

		lines[0] = Line(a1[0], a1[1], b1[0], b1[1]);
		lines[1] = Line(a2[0], a2[1], b2[0], b2[1]);
	}
	else
	{
		q::dvec apex = intersection(a1, a2, f, f + r1);

		if(side == KNAPSACK2_RIGHT)
		{
			lines[1] = Line(a1[0], a1[1], a2[0], a2[1]);
			lines[0] = Line(o[0], o[1], apex[0], apex[1]);
		}
		else
		{
			lines[0] = Line(a1[0], a1[1], a2[0], a2[1]);
			lines[1] = Line(o[0], o[1], apex[0], apex[1]);
		}
	}

	if(debug)
		cout << endl << lines[0] << "  " << lines[1] << endl;

	cur_facet++;

	IntersectionCut ic(Point(f[0], f[1]), 2, lines);
	cut->pi.clear();
	cut->pi_zero = -1;
	cut->pi.resize(row.c.pi.size());
	cut->depth = cur_facet-1;

	Point ray;
	rational alpha, rx, ry;

	for (int l = 0; l < row.c.pi.nz(); l++)
	{
		int j = row.c.pi.index(l);
		if (j == row.basic_var_index)
			continue;

		rx = -row.c.pi.value(l);
		ry = (l == r1_offset ? 1 : 0);

		if (row.is_integer[j] && l != r1_offset)
		{
			#ifdef INTERSECTION_CUT_USE_DOUBLE
				cut->pi.push(j, -ic.get_trivial_lifting_coefficient_double(
						rx.get_double(), ry.get_double()));
			#else
				cut->pi.push(j, -ic.get_trivial_lifting_coefficient(rx, ry));
			#endif
		}
		else
		{
			cut->pi.push(j, -ic.get_continuous_coefficient(rx, ry));
		}

		if(debug)
			cout << "r" << j << ": (" << rx << ", " << ry << ")  " << cut->pi[j] << endl;
	}
}
