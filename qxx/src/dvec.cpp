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
#include "qxx/debug.hpp"
#include "qxx/rational.hpp"
#include "qxx/dvec.hpp"
#include "qxx/svec.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {


// **************************************************************************
// 
// **************************************************************************
dvec::dvec()
{
}

dvec::dvec(const dvec &b)
	: d(b.d)
{
}

dvec::dvec(int n)
	: d(n)
{
}

dvec::dvec(const svec &b)
{
	b.spread(*this);
}


// **************************************************************************
// 
// **************************************************************************
int dvec::size() const
{
	return(d.size());
}

void dvec::resize(int n)
{
	d.resize(n);
}

void dvec::clear()
{
	d.clear();
}

mpq &dvec::operator[] (int i)
{
	Q_RANGE_CHECK(i, 0, size() - 1);
	return(d[i]);
}

const mpq &dvec::operator[] (int i) const
{
	Q_RANGE_CHECK(i, 0, size() - 1);
	return(d[i]);
}

// **************************************************************************
// 
// **************************************************************************
void dvec::set(const mpq &v)
{
	set(0, size() - 1, v);
}

void dvec::set(int i0, int i1, const mpq &v)
{
	for (int i = i0; i <= i1; i++)
		d[i] = v;
}

void dvec::assign(int n, const mpq &v)
{
	resize(n);
	set(v);
}

// **************************************************************************
// 
// **************************************************************************
dvec dvec::operator+(const dvec &b) const
{
	Q_MATCH(size(), b.size());

	int n = size();
	
	dvec r(n);
	
	for (int i = 0; i < n; i++)
		r.d[i] = d[i] + b.d[i];
	
	return(r);
}

dvec dvec::operator-(const dvec &b) const
{
	Q_MATCH(size(), b.size());
	
	int n = size();
	
	dvec r(n);
	
	for (int i = 0; i < n; i++)
		r.d[i] = d[i] - b.d[i];
	
	return(r);
}

mpq dvec::operator*(const dvec &b) const
{
	Q_MATCH(size(), b.size());
	
	int n = size();

	mpq r;
	
	for (int i = 0; i < n; i++)
		r += d[i] * b.d[i];
	
	return(r);
}

mpq dvec::operator*(const svec &b) const
{
	return(b.operator*(*this));
}


// **************************************************************************
// 
// **************************************************************************
const dvec &dvec::operator+() const
{
	return(*this);
}

dvec dvec::operator-() const
{
	int n = size();
	dvec r(n);
	
	for (int i = 0; i < n; i++)
		r.d[i] = -d[i];
	
	return(r);
}

dvec dvec::operator+(const mpq &v) const
{
	int n = size();
	
	dvec r(n);
	
	for (int i = 0; i < n; i++)
		r.d[i] = d[i] + v;
	
	return(r);
}

dvec dvec::operator-(const mpq &v) const
{
	int n = size();
	dvec r(n);
	
	for (int i = 0; i < n; i++)
		r.d[i] = d[i] - v;
	
	return(r);
}

dvec dvec::operator*(const mpq &v) const
{
	int n = size();
	dvec r(n);
	
	for (int i = 0; i < n; i++)
		r.d[i] = d[i] * v;
	
	return(r);
}

dvec dvec::operator/(const mpq &v) const
{
	int n = size();
	dvec r(n);
	
	for (int i = 0; i < n; i++)
		r.d[i] = d[i] / v;
	
	return(r);
}


// **************************************************************************
// 
// **************************************************************************
bool dvec::operator==(const dvec &b) const
{
	Q_MATCH(size(), b.size());
	
	int n = size();
	
	for (int i = 0; i < n; i++) {
		if (d[i] != b.d[i])
			return(false);
	}
	
	return(true);
}

bool dvec::operator!=(const dvec &b) const
{
	return(!((*this) == b));
}

// **************************************************************************
// 
// **************************************************************************
dvec &dvec::operator=(const dvec &b)
{
	d = b.d;
	return(*this);
}

dvec &dvec::operator=(const svec &b)
{
	b.spread(*this);
	return(*this);
}


// **************************************************************************
// 
// **************************************************************************
dvec &dvec::operator+=(const dvec &b)
{
	Q_MATCH(size(), b.size());
	
	int n = size();
	
	for (int i = 0; i < n; i++)
		d[i] += b.d[i];
	
	return(*this);
}

dvec &dvec::operator-=(const dvec &b)
{
	Q_MATCH(size(), b.size());

	int n = size();
	
	for (int i = 0; i < n; i++)
		d[i] -= b.d[i];
	
	return(*this);
}

dvec &dvec::operator+=(const mpq &v)
{
	int n = size();
	
	for (int i = 0; i < n; i++)
		d[i] += v;
	
	return(*this);
}

dvec &dvec::operator-=(const mpq &v)
{
	int n = size();
	
	for (int i = 0; i < n; i++)
		d[i] -= v;
	
	return(*this);
}

dvec &dvec::operator*=(const mpq &v)
{
	int n = size();
	
	for (int i = 0; i < n; i++)
		d[i] *= v;
	
	return(*this);
}

dvec &dvec::operator/=(const mpq &v)
{
	int n = size();
	
	for (int i = 0; i < n; i++)
		d[i] /= v;
	
	return(*this);
}


// **************************************************************************
// 
// **************************************************************************
void dvec::fdump(FILE *f, const std::string &name) const
{
	int n = size();
	
	if (name.length())
		gmp_fprintf(f, "%s = [", name.c_str());
	
	for (int i = 0; i < n; i++) {
		d[i].fdump(f);
	}
	
	if (name.length())
		gmp_fprintf(f, " ]';\n");
	else
		gmp_fprintf(f, ";\n");
}

void dvec::dump(const std::string &name) const
{
	fdump(stdout, name);
}

// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const dvec &v)
{
	int n = v.size();
	
	for (int i = 0; i < n; i++)
		os << " " << v[i];
	
	return(os);
}

// **************************************************************************
// 
// **************************************************************************

} // namespace

