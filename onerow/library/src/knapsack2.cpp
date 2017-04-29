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

#include <cstdio>
#include <cstdlib>

#include <qxx/rational.hpp>
#include <qxx/dlu.hpp>

#include <onerow/knapsack2.hpp>
#include <onerow/geometry.hpp>

using std::endl;
using std::cout;


// **************************************************************************
// 
// **************************************************************************
//#define DEBUG(s...)	gmp_printf(s)
#define DEBUG(s...)

// **************************************************************************
// 
// **************************************************************************
Knapsack2::Knapsack2()
{
}

Knapsack2::Knapsack2(rational f, rational r1)
{
	eval(f, r1);
}

Knapsack2::~Knapsack2()
{
}

// **************************************************************************
// 
// **************************************************************************
void Knapsack2::clear()
{
	list.clear();
}


// **************************************************************************
// 
// **************************************************************************
void Knapsack2::push(int side,
	const q::dvec &l, const q::dvec &u, const q::dvec &o)
{
	Knapsack2Vertex v;
	
	v.side = side;
	v.lower = l;
	v.upper = u;
	v.opposed = o;
	v.lower.resize(2);
	v.upper.resize(2);
	v.opposed.resize(2);
	
	list.push_back(v);
}

//int hit_count = 0;
//int miss_count = 0;

// **************************************************************************
// 
// **************************************************************************
void Knapsack2::eval(rational fx, rational r1x)
{
	clear();
	
	rational fr1frac, fr1floor, hslope;
	q::dvec a(3), b(3), p(3), f(3), g(3), h(3), r1(3), wy(3), wa(3), wb(3);
	q::dmat w(3, 3), u(3, 3), ta(3, 3), tb(3, 3), tx(3, 3);
	
	a[0] = 0; a[1] = 0; a[2] = 1;
	b[0] = 1; b[1] = 0; b[2] = 1;
	p[2] = 1;
	f[0] = fx.frac(); f[1] = 0; f[2] = 1;
	g[2] = 1;
	h[2] = 1;
	r1[0] = r1x; r1[1] = 1; r1[2] = 0;
	
	w.set_identity();
	w(0, 2) = fx.floor();
	tb(0, 1) = 1; tb(0, 2) = 1; tb(1, 2) = 1;
	tb(2, 0) = tb(2, 1) = tb(2, 2) = 1;
	ta(0, 1) = 1; ta(1, 2) = 1;
	ta(2, 0) = ta(2, 1) = ta(2, 2) = 1;
	
	int it = 0;
	
	wa = w * a;
	DEBUG("-: A vertex (0, 0)\t\t--> v A  (%Qd, %Qd)\n",
		wa[0].v, wa[1].v);
	wb = w * b;
	DEBUG("-: B vertex (1, 0)\t\t--> v  B (%Qd, %Qd)\n",
		wb[0].v, wb[1].v);

	while (1) {
		DEBUG("%d: -----------------------------\n", it);
		DEBUG("%d: f = (%Qd, %Qd, %Qd)\n",
			it, f[0].v, f[1].v, f[2].v);
		DEBUG("%d: r1 = (%Qd, %Qd, %Qd)\n",
			it, r1[0].v, r1[1].v, r1[2].v);
		DEBUG("%d: a = (%Qd, %Qd, %Qd)\n",
			it, a[0].v, a[1].v, a[2].v);
		DEBUG("%d: b = (%Qd, %Qd, %Qd)\n",
			it, b[0].v, b[1].v, b[2].v);

		// Step 1
		fr1frac = (f[0] + r1[0]).frac();
		fr1floor = (f[0] + r1[0]).floor();
		
		if (f[0] == fr1frac) {
			wy = w * r1;
			DEBUG("%d: AB ray (%Qd, %Qd)"
				"\t\t--> r AB (%Qd, %Qd)\n",
				it, r1[0].v, r1[1].v, wy[0].v, wy[1].v);
			push(KNAPSACK2_RAY, wa, wy, wb);
			break;
		} else if (f[0] < fr1frac) {
			DEBUG("%d: hit right\n", it);

			g[1] = (q::mpq(1) - f[0]) / (fr1frac - f[0]);
			g[0] = q::mpq(1) + fr1floor * g[1];
			
			b[1] = g[1].floor();
			b[0] = q::mpq(1) + fr1floor * b[1];
			
			p[1] = b[1] + 1;
			p[0] = q::mpq(1) + fr1floor * p[1];
			
			if (b[1] == g[1]) {
				wy = w * b;
				DEBUG("%d: AB vertex (%Qd, %Qd)"
					"\t\t--> v AB (%Qd, %Qd)\n",
					it, b[0].v, b[1].v, wy[0].v, wy[1].v);
				push(KNAPSACK2_BOTH, wa, wy, wb);
				break;
			}
			
			tx(0, 0) = 0;
			tx(1, 0) = 0;
			tx(2, 0) = 1;
			tx.set_col(1, b);
			tx.set_col(2, p);
			
			u = tb * tx.inv();
			
			//tb.dump("tb");
			//tx.dump("tx");
			//u.dump("u");

			wy = w * b;

			push(KNAPSACK2_RIGHT, wb, wy, wa);
			
			wb = wy;
			
			DEBUG("%d: B vertex (%Qd, %Qd)"
				"\t\t--> v  B (%Qd, %Qd)"
				"   opposed (%Qd, %Qd)\n",
				it, b[0].v, b[1].v, wb[0].v, wb[1].v,
				wa[0].v, wa[1].v);
			
			hslope = b[0] / b[1];
			h[1] = f[0] / (hslope - r1[0]);
			h[0] = hslope * h[1];
		} else {
			DEBUG("%d: hit left\n", it);

			g[1] = f[0] / (f[0] - fr1frac);
			g[0] = fr1floor * g[1];
			
			a[1] = g[1].floor();
			a[0] = fr1floor * a[1];
			
			p[1] = a[1] + 1;
			p[0] = fr1floor * p[1];
			
			if (a[1] == g[1]) {
				wy = w * a;
				DEBUG("%d: AB vertex (%Qd, %Qd)"
					"\t\t--> v AB (%Qd, %Qd)\n",
					it, a[0].v, a[1].v, wy[0].v, wy[1].v);
				push(KNAPSACK2_BOTH, wa, wy, wb);
				break;
			}
			
			tx.set_col(0, a);
			tx(0, 1) = 1;
			tx(1, 1) = 0;
			tx(2, 1) = 1;
			tx.set_col(2, p);
			
			u = ta * tx.inv();

			//tb.dump("ta");
			//u.dump("u");

			wy = w * a;

			push(KNAPSACK2_LEFT, wa, wy, wb);

			wa = wy;
			
			DEBUG("%d: A vertex (%Qd, %Qd)"
				"\t\t--> v A  (%Qd, %Qd)"
				"   opposed (%Qd, %Qd)\n",
				it, a[0].v, a[1].v, wa[0].v, wa[1].v,
				wb[0].v, wb[1].v);
			
			hslope = (a[0] - 1) / a[1];
			h[1] = (f[0] - 1) / (hslope - r1[0]);
			h[0] = q::mpq(1) + (hslope * h[1]);
		}
		
		a = u * a;
		b = u * b;
		f = u * h;
		r1 = u * r1;
		r1[0] /= r1[1];
		r1[1] = 1;
		
		w = w * u.inv();
		
		it++;
	}
}

