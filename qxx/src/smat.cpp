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
#include <cstdio>
#include "qxx/debug.hpp"
#include "qxx/rational.hpp"
#include "qxx/dvec.hpp"
#include "qxx/slu.hpp"
#include "qxx/smat.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {


// **************************************************************************
// 
// **************************************************************************
smat::smat()
	: m(0), n(0)
{
}

smat::smat(const smat &b)
	: m(b.m), n(b.n), d(b.d)
{
}

smat::smat(int m0, int n0)
	: m(m0), n(n0), d(m0)
{
	for (int i = 0; i < m; i++)
		d[i].resize(n);
}


int smat::rows() const
{
	return(m);
}

int smat::cols() const
{
	return(n);
}

void smat::resize(int m0, int n0)
{
	m = m0;
	n = n0;
	
	d.resize(m);
	for (int i = 0; i < m; i++)
		d[i].resize(n);
}

void smat::clear()
{
	m = n = 0;
	d.clear();
}

// **************************************************************************
// 
// **************************************************************************
const mpq &smat::get(int i, int j) const
{
	Q_RANGE_CHECK(i, 0, m - 1);
	Q_RANGE_CHECK(j, 0, n - 1);

	return(d[i][j]);
}

void smat::set(int i, int j, const mpq &v)
{
	Q_RANGE_CHECK(i, 0, m - 1);
	Q_RANGE_CHECK(j, 0, n - 1);

	d[i][j] = v;
}

void smat::set_zero()
{
	for (int i = 0; i < m; i++)
		d[i].set_zero();
}

void smat::set_identity()
{
	mpq one(1);
	
	for (int i = 0; i < m; i++) {
		d[i].set_zero();
		d[i][i] = one;
	}
}


// **************************************************************************
// 
// **************************************************************************
void smat::gather(const dmat &src)
{
	resize(src.rows(), src.cols());
	
	for (int i = 0; i < m; i++)
		d[i].gather(src[i]);
}

void smat::spread(dmat &r) const
{
	r.resize(m, n);
	
	for (int i = 0; i < m; i++)
		d[i].spread(r[i]);
}

dmat smat::dense() const
{
	dmat r;
	
	spread(r);
	
	return(r);
}

// **************************************************************************
// 
// **************************************************************************
svec smat::get_col(int j) const
{
	svec r(m);
	
	for (int i = 0; i < m; i++)
		r.push(i, d[i][j]);
	
	return(r);
}

void smat::set_col(int j, const svec &v)
{
	Q_MATCH(v.size(), m);
	
	for (int i = 0; i < m; i++)
		d[i][j] = v[i];
}

const svec &smat::get_row(int i) const
{
	return(d[i]);
}

void smat::set_row(int i, const svec &v)
{
	Q_MATCH(v.size(), n);
	
	d[i] = v;
}


// **************************************************************************
// 
// **************************************************************************
svec &smat::operator[] (int i)
{
	return(d[i]);
}

const svec &smat::operator[] (int i) const
{
	return(d[i]);
}

// **************************************************************************
// 
// **************************************************************************
svec_ref smat::operator() (int i, int j)
{
	Q_RANGE_CHECK(i, 0, m - 1);
	Q_RANGE_CHECK(j, 0, n - 1);

	return(d[i][j]);
}

const mpq &smat::operator() (int i, int j) const
{
	Q_RANGE_CHECK(i, 0, m - 1);
	Q_RANGE_CHECK(j, 0, n - 1);

	return(d[i][j]);
}


// **************************************************************************
// 
// **************************************************************************
smat &smat::operator=(const smat &a)
{
	m = a.m;
	n = a.n;
	d = a.d;
	
	return(*this);
}

// **************************************************************************
// 
// **************************************************************************
dvec smat::operator*(const dvec &b) const
{
	Q_MATCH(n, b.size());
	
	dvec x(m);
	mpq v;
	
	for (int i = 0; i < m; i++) {
		v = 0;
		
		const svec &r = d[i];
		for (int l = 0; l < r.nz(); l++)
			v += r.value(l) * b[r.index(l)];
		
		x[i] = v;
	}
	
	return(x);
}

// **************************************************************************
// 
// **************************************************************************
void smat::set_transpose(const smat &a)
{
	resize(a.n, a.m);
	
	for (int i = 0; i < a.n; i++) {
		const svec &ar = a[i];
		
		for (int l = 0; l < ar.nz(); l++)
			d[ar.index(l)].push(i, ar.value(l));
	}
}

// **************************************************************************
// 
// **************************************************************************
smat smat::t() const
{
	smat a;
	
	a.set_transpose(*this);
	
	return(a);
}


// **************************************************************************
// 
// **************************************************************************
smat smat::inv() const
{
	Q_MATCH(m, n);
	
	slu lu(*this);
	smat g;
	dvec e;
	svec x;
	
	g.resize(n, n);
	e.assign(n, 0);
	
	for (int i = 0; i < n; i++) {
		e[i] = 1;

		x = lu.solve_xA(e);
		
		g.set_row(i, x);
		
		e[i] = 0;
	}
	
	return(g);
}


// **************************************************************************
// 
// **************************************************************************
void smat::fdump(FILE *f, const std::string &name) const
{
	if (name.length())
		gmp_fprintf(f, "%s = [\n", name.c_str());
	else
		gmp_fprintf(f, "[\n");

	for (int i = 0; i < m; i++)
		d[i].fdump(f);

	gmp_fprintf(f, "];\n");
}

void smat::dump(const std::string &name) const
{
	fdump(stdout, name);
}

// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const smat &a)
{
	for (int i = 0; i < a.rows(); i++)
		os << a[i] << std::endl;
	
	return(os);
}

// **************************************************************************
// 
// **************************************************************************
}




