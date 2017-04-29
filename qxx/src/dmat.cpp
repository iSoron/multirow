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
#include "qxx/dlu.hpp"
#include "qxx/dmat.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {


// **************************************************************************
// 
// **************************************************************************
dmat::dmat()
	: m(0), n(0)
{
}

dmat::dmat(const dmat &b)
	: m(b.m), n(b.n), d(b.d)
{
}

dmat::dmat(int m0, int n0)
	: m(m0), n(n0), d(m0)
{
	for (int i = 0; i < m; i++)
		d[i].resize(n);
}


int dmat::rows() const
{
	return(m);
}

int dmat::cols() const
{
	return(n);
}

void dmat::resize(int m0, int n0)
{
	m = m0;
	n = n0;
	
	d.resize(m);
	for (int i = 0; i < m; i++)
		d[i].resize(n);
}

void dmat::clear()
{
	m = n = 0;
	d.clear();
}

// **************************************************************************
// 
// **************************************************************************
const mpq &dmat::get(int i, int j) const
{
	Q_RANGE_CHECK(i, 0, m - 1);
	Q_RANGE_CHECK(j, 0, n - 1);
	
	return(d[i][j]);
}

void dmat::set(int i, int j, const mpq &v)
{
	Q_RANGE_CHECK(i, 0, m - 1);
	Q_RANGE_CHECK(j, 0, n - 1);

	d[i][j] = v;
}

void dmat::set(const mpq &v)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			d[i][j] = v;
}

void dmat::set_identity()
{
	mpq zero(0), one(1);
	
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			d[i][j] = (i == j) ? one : zero;
}


// **************************************************************************
// 
// **************************************************************************
dvec dmat::get_col(int j) const
{
	Q_RANGE_CHECK(j, 0, n - 1);

	dvec r(m);
	
	for (int i = 0; i < m; i++)
		r[i] = d[i][j];
	
	return(r);
}

void dmat::set_col(int j, const dvec &v)
{
	Q_RANGE_CHECK(j, 0, n - 1);
	Q_MATCH(v.size(), m);
	
	for (int i = 0; i < m; i++)
		d[i][j] = v[i];
}

const dvec &dmat::get_row(int i) const
{
	Q_RANGE_CHECK(i, 0, m - 1);

	return(d[i]);
}

void dmat::set_row(int i, const dvec &v)
{
	Q_RANGE_CHECK(i, 0, m - 1);
	Q_MATCH(v.size(), n);
	
	d[i] = v;
}


// **************************************************************************
// 
// **************************************************************************
dvec &dmat::operator[] (int i)
{
	Q_RANGE_CHECK(i, 0, m - 1);

	return(d[i]);
}

const dvec &dmat::operator[] (int i) const
{
	Q_RANGE_CHECK(i, 0, m - 1);

	return(d[i]);
}

mpq &dmat::operator() (int i, int j)
{
	Q_RANGE_CHECK(i, 0, m - 1);
	Q_RANGE_CHECK(j, 0, n - 1);

	return(d[i][j]);
}

const mpq &dmat::operator() (int i, int j) const
{
	Q_RANGE_CHECK(i, 0, m - 1);
	Q_RANGE_CHECK(j, 0, n - 1);

	return(d[i][j]);
}


// **************************************************************************
// 
// **************************************************************************
const dmat &dmat::operator+() const
{
	return(*this);
}

dmat dmat::operator-() const
{
	dmat r(m, n);
	
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			r.d[i][j] = -d[i][j];
	
	return(r);
}

dmat dmat::operator+(const dmat &b) const
{
	Q_MATCH(m, b.m);
	Q_MATCH(n, b.n);
	
	dmat r(m, n);
	
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			r.d[i][j] = d[i][j] + b.d[i][j];
	
	return(r);
}

dmat dmat::operator-(const dmat &b) const
{
	Q_MATCH(m, b.m);
	Q_MATCH(n, b.n);
	
	dmat r(m, n);
	
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			r.d[i][j] = d[i][j] - b.d[i][j];
	
	return(r);
}

dmat dmat::operator*(const dmat &b) const
{
	Q_MATCH(n, b.m);
	
	dmat r(m, b.n);
	mpq v;
	
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < b.n; j++) {
			v = 0;
			
			for (int k = 0; k < n; k++)
				v += d[i][k] * b.d[k][j];
			
			r.d[i][j] = v;
		}
	}
	
	return(r);
}

dvec dmat::operator*(const dvec &b) const
{
	Q_MATCH(n, b.size());
	
	dvec r(m);
	mpq v;
	
	for (int i = 0; i < m; i++) {
		v = 0;
		for (int j = 0; j < n; j++)
			v += d[i][j] * b[j];
		r[i] = v;
	}
	
	return(r);
}

dmat dmat::operator*(const mpq &v) const
{
	dmat r(m, n);
	
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			r.d[i][j] = d[i][j] * v;
	
	return(r);
}

dmat dmat::operator/(const mpq &v) const
{
	dmat r(m, n);
	
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			r.d[i][j] = d[i][j] / v;
	
	return(r);
}


// **************************************************************************
// 
// **************************************************************************
dmat dmat::inv() const
{
	Q_MATCH(m, n);
	
	dvec e(n);
	dmat r(n, n);
	dlu lu(*this);
	
	for (int i = 0; i < n; i++) {
		e[i] = 1;
		
		r.set_col(i, lu.solve_Ax(e));
		
		//e.dump(stdout, "e");
		//r.get_col(i).dump(stdout, "x");
		
		e[i] = 0;
	}
	
	return(r);
}

dmat dmat::t() const
{
	dmat r(m, n);
	
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			r.d[i][j] = d[j][i];
	
	return(r);
}

mpq dmat::det() const
{
	Q_MATCH(m, n);
	
	if (n == 1)
		return(d[0][0]);
	
	if (n == 2)
		return(d[0][0] * d[1][1] - d[1][0] * d[0][1]);
	
	if (n == 3)
		return(	  d[0][0] * d[1][1] * d[2][2]
			+ d[0][1] * d[1][2] * d[2][0]
			+ d[0][2] * d[1][0] * d[2][1]
			- d[2][0] * d[1][1] * d[0][2]
			- d[2][1] * d[1][2] * d[0][0]
			- d[2][2] * d[1][0] * d[0][1] );
	
	dlu lu(*this);
	return(lu.det());
}

// **************************************************************************
// 
// **************************************************************************
bool dmat::operator==(const dmat &b) const
{
	Q_MATCH(m, b.m);
	Q_MATCH(n, b.n);
	
	for (int i = 0; i < m; i++) {
		if (d[i] != b.d[i])
			return(false);
	}
	
	return(true);
}

bool dmat::operator!=(const dmat &b) const
{
	return(!((*this) == b));
}

// **************************************************************************
// 
// **************************************************************************
dmat &dmat::operator=(const dmat &b)
{
	m = b.m;
	n = b.n;
	d = b.d;
	return(*this);
}

dmat &dmat::operator=(const dvec &v)
{
	resize(v.size(), 1);
	
	for (int i = 0; i < m; i++)
		d[i][0] = v[i];
	
	return(*this);
}


dmat &dmat::operator+=(const dmat &b)
{
	Q_MATCH(m, b.m);
	Q_MATCH(n, b.n);
	
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			d[i][j] += b.d[i][j];
	
	return(*this);
}

dmat &dmat::operator-=(const dmat &b)
{
	Q_MATCH(m, b.m);
	Q_MATCH(n, b.n);
	
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			d[i][j] -= b.d[i][j];
	
	return(*this);
}

dmat &dmat::operator*=(const dmat &b)
{
	dmat r;
	
	r = (*this) * b;
	(*this) = r;
	return(*this);
}

dmat &dmat::operator*=(const mpq &v)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			d[i][j] *= v;
	
	return(*this);
}

dmat &dmat::operator/=(const mpq &v)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			d[i][j] /= v;
	
	return(*this);
}


// **************************************************************************
// 
// **************************************************************************
void dmat::fdump(FILE *f, const std::string &name) const
{
	if (name.length())
		gmp_fprintf(f, "%s = [\n", name.c_str());
	else
		gmp_fprintf(f, "[\n");
	
	for (int i = 0; i < m; i++)
		d[i].fdump(f);
	fprintf(f, "];\n");
}

void dmat::dump(const std::string &name) const
{
	fdump(stdout, name);
}

// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const dmat &a)
{
	for (int i = 0; i < a.rows(); i++)
		os << a[i] << std::endl;
	
	return(os);
}

// **************************************************************************
// 
// **************************************************************************
}




