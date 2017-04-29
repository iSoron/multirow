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
#include "qxx/strprintf.hpp"
#include "qxx/rational.hpp"
#include "qxx/slu.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {


// **************************************************************************
// 
// **************************************************************************
perm::perm()
{
}

perm::~perm()
{
}

int perm::size() const
{
	return(fwd.size());
}

void perm::resize(int n)
{
	fwd.resize(n);
	bwd.resize(n);
}

void perm::clear()
{
	fwd.clear();
	bwd.clear();
}

void perm::id()
{
	int n = fwd.size();
	
	for (int i = 0; i < n; i++)
		fwd[i] = i;
	
	for (int k = 0; k < n; k++)
		bwd[k] = k;
}

void perm::id(int n)
{
	resize(n);
	id();
}

void perm::pivot(int k, int i)
{
	int t = fwd[i];
	int u = bwd[k];
	
	fwd[i] = k;
	fwd[u] = t;
	bwd[k] = i;
	bwd[t] = u;
}


// **************************************************************************
// 
// **************************************************************************
slu::slu(const smat &a)
{
	factorize(a);
}

slu::~slu()
{
}

// **************************************************************************
// 
// **************************************************************************
void slu::pivot(int p)
{
	double best_v = n * n;
	int best_i = -1;
	int best_j = -1;
	
	for (int k = p; k < n; k++) {
		int i = row.bwd[k];
		int inz = w[i].nz();
		
		for (int l = 0; l < inz; l++) {
			int j = w[i].index(l);
			int jnz = cnz[j];
			int kj = col.fwd[j];
			double v = (double)(inz - 1) * (jnz - 1);
			
			if ((kj >= p) && (v < best_v)) {
				best_v = v;
				best_i = i;
				best_j = j;
			}
		}
	}
	
	assert((best_i >= 0) && "Singular matrix");
	
	//best_i = best_j = p;
	row.pivot(p, best_i);
	col.pivot(p, best_j);
}

// **************************************************************************
// 
// **************************************************************************
void slu::factorize(const smat &a)
{
	// check
	assert(a.rows() == a.cols());
	
	// init
	n = a.rows();
	w = a;
	lt.resize(n, n);
	u.resize(n, n);
	row.id(n);
	col.id(n);
	
	cnz.assign(n, 0);

	for (int i = 0; i < n; i++) {
		for (int l = 0; l < w[i].nz(); l++)
			cnz[w[i].index(l)]++;
	}
	
	// main loop
	svec pr;
	mpq pv, f;
	std::vector<bool> hit;
	
	for (int p = 0; p < n; p++) {
		//w.dense().dump(strprintf("w%d", p).c_str());
	
		pivot(p);
		
		int pi = row.bwd[p];
		int pj = col.bwd[p];
		
		// get pivot row
		pr = w[pi];
		pr.offs_build();
		pv = pr[pj];
		for (int pl = 0; pl < pr.nz(); pl++)
			cnz[pr.index(pl)]--;
		
		// update u, lt, w
		u[p] = w[pi];
		lt[p] = w.get_col(pj) / pv;
		w[pi].set_zero();
		
		// elimination
		for (int k = p + 1; k < n; k++) {
			int i = row.bwd[k];
			svec &ir = w[i];
			
			f = -ir[pj] / pv;
			
			if (f.sign() == 0)
				continue;
			
			hit.assign(pr.nz(), false);

			for (int l = 0; l < ir.nz(); l++) {
				int j = ir.index(l);
				int pl = pr.locate(j);
				
				if (pl == -1)
					continue;
				
				hit[pl] = true;

				ir.value(l) += f * pr.value(pl);

				if (ir.value(l).sign() == 0) {
					cnz[j]--;
					ir.remove(l);
					l--;
				}
			}
			
			for (int pl = 0; pl < pr.nz(); pl++) {
				if (hit[pl])
					continue;
				int j = pr.index(pl);
				
				ir.push(j, f * pr.value(pl));
				cnz[j]++;
			}
		}
	}
	
	// postprocessing
	for (int p = 0; p < n; p++) {
		svec &pr = u[p];
		
		for (int l = 0; l < pr.nz(); l++)
			pr.index(l) = col.fwd[pr.index(l)];
	}


	for (int p = 0; p < n; p++) {
		svec &pc = lt[p];
		
		for (int l = 0; l < pc.nz(); l++)
			pc.index(l) = row.fwd[pc.index(l)];
	}
	
	l.set_transpose(lt);
	ut.set_transpose(u);
}


// **************************************************************************
// 
// **************************************************************************
dvec slu::solve_gen(const smat &l, const smat &u,
	const perm &row, const perm &col, const dvec &b) const
{
	dvec x(n), y(n);
	mpq v, diag;
	
	// Ly = b
	for (int k = 0; k < n; k++) {
		v = b[row.bwd[k]];
		
		const svec &lr = l[k];
		
		for (int l = 0; l < lr.nz(); l++) {
			int j = lr.index(l);
			
			if (j == k)
				diag = lr.value(l);
			else
				v -= lr.value(l) * y[j];
		}
		
		y[k] = v / diag;
	}
	
	// Ux = y
	for (int k = n - 1; k != -1; k--) {
		v = y[k];
		
		const svec &ur = u[k];
		
		for (int l = 0; l < ur.nz(); l++) {
			int j = ur.index(l);
			
			if (j == k)
				diag = ur.value(l);
			else
				v -= ur.value(l) * y[j];
		}
		
		y[k] = v / diag;
	}
	
	// postprocessing
	for (int j = 0; j < n; j++)
		x[j] = y[col.fwd[j]];
	
	return(x);
}

dvec slu::solve_Ax(const dvec &b) const
{
	return(solve_gen(l, u, row, col, b));
}

dvec slu::solve_xA(const dvec &b) const
{
	return(solve_gen(ut, lt, col, row, b));
}

// **************************************************************************
// 
// **************************************************************************
};

