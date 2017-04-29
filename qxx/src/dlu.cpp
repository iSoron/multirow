/*
    This file is part of qxx -- matrix algebra in exact arithmetic
    Copyright (C) 2013-2014  Laurent Poirrier

    libp is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pxx.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <cassert>
#include "qxx/dlu.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {


// **************************************************************************
// 
// **************************************************************************
dlu::dlu(const dmat &a)
	: n(0)
{
	factorize(a);
}

dlu::~dlu()
{
}

// **************************************************************************
// 
// **************************************************************************
void dlu::pivot(int p)
{
	int best_v = -1;
	int best_k = -1;
	
	//printf("pivot(%d):\n", p);
	
	for (int k = p; k < n; k++) {
		int i = row_bwd[k];
		
		//gmp_printf("\tchoice [%d] row %d: %Qd\n",
		//	k, i, g(i, p).v);
		
		if (g(i, p).sign()) {
			int v = g(i, p).bits();
			
			if ((best_v < 0) || (v < best_v)) {
				best_k = k;
				best_v = v;
			}
		}
	}
	
	assert((best_k >= 0) && "Singular matrix");
	
	int i = row_bwd[best_k];
	row_bwd[best_k] = row_bwd[p];
	row_bwd[p] = i;
}

// **************************************************************************
// 
// **************************************************************************
void dlu::factorize(const dmat &a)
{
	// check
	assert(a.rows() == a.cols());
	
	// init
	n = a.rows();
	g = a;
	row_bwd.resize(n);
	for (int p = 0; p < n; p++)
		row_bwd[p] = p;
	
	// main loop
	mpq pv, f;
	
	for (int p = 0; p < n; p++) {
		pivot(p);

		int pi = row_bwd[p];
		pv = g(pi, p);
		
		for (int k = p + 1; k < n; k++) {
			int i = row_bwd[k];
			f = g(i, p) / pv;
			
			for (int j = p + 1; j < n; j++)
				g(i, j) -= f * g(pi, j);
		}
		
		//printf("g%d = ", p);
		//g.dump(stdout);
	}

	//printf("row_bwd = [");
	//for (int p = 0; p < n; p++)
	//	printf(" %d", row_bwd[p]);
	//printf(" ];\n");
}		

// **************************************************************************
// 
// **************************************************************************
dmat dlu::L() const
{
	dmat l(n, n);
	
	for (int k = 0; k < n; k++) {
		for (int j = 0; j < k; j++)
			l(k, j) = g(row_bwd[k], j) / g(row_bwd[j], j);
		l(k, k) = 1;
	}
	
	return(l);
}

dmat dlu::U() const
{
	dmat u(n, n);
	
	for (int k = 0; k < n; k++) {
		for (int j = k; j < n; j++)
			u(k, j) = g(row_bwd[k], j);
	}
	
	return(u);
}

// **************************************************************************
// 
// **************************************************************************
dvec dlu::solve_Ax(const dvec &b) const
{
	assert(b.size() == n);
	
	dvec x(n);
	mpq v;
	
	int i0 = row_bwd[0];
	x[0] = b[i0] / g(i0, 0);
	
	for (int k = 1; k < n; k++) {
		int i = row_bwd[k];
		v = 0;
		
		for (int j = 0; j < k; j++)
			v += g(i, j) * x[j];
		
		x[k] = (b[i] - v) / g(i, k);
	}
	
	for (int k = n - 2; k != -1; k--) {
		int i = row_bwd[k];
		v = 0;
		
		for (int j = k + 1; j < n; j++)
			v += g(i, j) * x[j];
		
		x[k] -= v / g(i, k);
	}
	
	return(x);
}

// **************************************************************************
// 
// **************************************************************************
dvec dlu::solve_xA(const dvec &b) const
{
	assert(b.size() == n);
	
	dvec x(n);
	mpq v;
	
	int i0 = row_bwd[0];
	x[i0] = b[0] / g(i0, 0);
	
	for (int k = 1; k < n; k++) {
		int i = row_bwd[k];
		v = 0;

		for (int l = 0; l < k; l++) {
			int il = row_bwd[l];
			v += g(il, k) * x[il];
		}
		
		x[i] = (b[k] - v) / g(i, k);
	}
	
	for (int k = n - 2; k != -1; k--) {
		int i = row_bwd[k];
		v = 0;
		
		for (int l = k + 1; l < n; l++) {
			int il = row_bwd[l];
			v += g(il, k) * x[il];
		}
		
		x[i] -= v / g(i, k);
	}
	
	return(x);
}

// **************************************************************************
// 
// **************************************************************************
mpq dlu::det() const
{
	mpq v;
	
	v = 1;
	for (int k = 0; k < n; k++) {
		int i = row_bwd[k];
		
		v *= g(i, k);
	}
	
	return(v);
}

// **************************************************************************
// 
// **************************************************************************
};


