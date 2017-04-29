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
#include "qxx/svec.hpp"
#include <algorithm>

// **************************************************************************
// 
// **************************************************************************
namespace q {


// **************************************************************************
// svec_ref: construction and destruction
// **************************************************************************
svec_ref::svec_ref(svec &vec, int i)
	: vector(vec), index(i)
{
}

svec_ref::svec_ref(const svec_ref &src)
	: vector(src.vector), index(src.index)
{
}

svec_ref::~svec_ref()
{
}


// **************************************************************************
// svec_ref: methods
// **************************************************************************
const mpq &svec_ref::get() const
{
	return(vector.get(index));
}

// **************************************************************************
// svec_ref: assignment operators
// **************************************************************************
svec_ref::operator mpq() const
{
	return(get());
}

const mpq &svec_ref::operator=(const svec_ref &b) const
{
	return(vector.set(index, b.get()));
}

const mpq &svec_ref::operator=(const mpq &b) const
{
	return(vector.set(index, b));
}

const mpq &svec_ref::operator+=(const mpq &b) const
{
	return(vector.set(index, vector.get(index) + b));
}

const mpq &svec_ref::operator-=(const mpq &b) const
{
	return(vector.set(index, vector.get(index) - b));
}

const mpq &svec_ref::operator*=(const mpq &b) const
{
	return(vector.set(index, vector.get(index) * b));
}

const mpq &svec_ref::operator/=(const mpq &b) const
{
	return(vector.set(index, vector.get(index) / b));
}

// **************************************************************************
// svec_ref: read-only operators
// **************************************************************************
const mpq &svec_ref::operator+() const
{
	return(get());
}

mpq svec_ref::operator-() const
{
	return(-get());
}

mpq svec_ref::operator+(const mpq &b) const
{
	return(get() + b);
}

mpq svec_ref::operator-(const mpq &b) const
{
	return(get() - b);
}

mpq svec_ref::operator*(const mpq &b) const
{
	return(get() * b);
}

mpq svec_ref::operator/(const mpq &b) const
{
	return(get() / b);
}


bool svec_ref::operator<(const mpq &b) const
{
	return(get() < b);
}

bool svec_ref::operator>(const mpq &b) const
{
	return(get() > b);
}

bool svec_ref::operator<=(const mpq &b) const
{
	return(get() <= b);
}

bool svec_ref::operator>=(const mpq &b) const
{
	return(get() >= b);
}

bool svec_ref::operator==(const mpq &b) const
{
	return(get() == b);
}

bool svec_ref::operator!=(const mpq &b) const
{
	return(get() != b);
}

// **************************************************************************
// svec_el
// **************************************************************************
svec_el::svec_el()
{
}

svec_el::svec_el(int idx, const mpq &v)
	: index(idx), value(v)
{
}

svec_el::svec_el(const svec_el &e)
	: index(e.index), value(e.value)
{
}

svec_el::~svec_el()
{
}

// **************************************************************************
// svec
// **************************************************************************
const mpq svec::zero;

svec::svec()
	: n(0)
{
}

svec::svec(const svec &src)
	: n(src.n), d(src.d)
{
}

svec::svec(int n0)
	: n(n0)
{
}

svec::svec(const dvec &src)
	: n(src.size())
{
	gather(src);
}

svec::~svec()
{
}


// **************************************************************************
// 
// **************************************************************************
void svec::clear()
{
	n = 0;
	d.clear();
}


int svec::size() const
{
	return(n);
}

void svec::resize(int n0)
{
	n = n0;
}

int svec::nz() const
{
	return(d.size());
}

// **************************************************************************
// 
// **************************************************************************
bool svec::offs_active() const
{
	return(offs.size() != 0);
}

void svec::offs_build()
{
	if (offs_active())
		return;
	
	offs.resize(n);
	
	for (int i = 0; i < n; i++)
		offs[i] = 0;
	
	for (int l = 0; l < nz(); l++)
		offs[d[l].index] = l + 1;
}

void svec::offs_clear()
{
	offs.resize(0);
}


// **************************************************************************
// 
// **************************************************************************
void svec::set_zero()
{
	d.clear();
	
	if (offs_active()) {
		for (int i = 0; i < n; i++)
			offs[i] = 0;
	}
}

// **************************************************************************
// 
// **************************************************************************
void svec::gather(const dvec &src)
{
	n = src.size();
	d.clear();
	for (int i = 0; i < n; i++)
		push(i, src[i]);
}

void svec::spread(dvec &r) const
{
	r.resize(n);
	r.set(zero);
	for (int l = 0; l < nz(); l++)
		r[d[l].index] = d[l].value;
}

dvec svec::dense() const
{
	dvec r;
	
	spread(r);
	return(r);
}

// **************************************************************************
// 
// **************************************************************************
int svec::locate(int idx) const
{
	Q_RANGE_CHECK(idx, 0, n - 1);
	
	if (offs_active())
		return(offs[idx] - 1);
	
	for (int l = 0; l < nz(); l++) {
		if (d[l].index == idx)
			return(l);
	}
	
	return(-1);
}

// **************************************************************************
// 
// **************************************************************************
void svec::remove(int l)
{
	Q_RANGE_CHECK(l, 0, nz() - 1);

	int k = nz() - 1;
	
	if (offs_active())
		offs[d[l].index] = 0;
	
	if (l < k)
		d[l] = d[k];
		
	d.resize(k);
}

// **************************************************************************
// 
// **************************************************************************
int &svec::index(int offset)
{
	Q_RANGE_CHECK(offset, 0, nz() - 1);

	return(d[offset].index);
}

mpq &svec::value(int offset)
{
	Q_RANGE_CHECK(offset, 0, nz() - 1);

	return(d[offset].value);
}


int svec::index(int offset) const
{
	Q_RANGE_CHECK(offset, 0, nz() - 1);

	return(d[offset].index);
}

const mpq &svec::value(int offset) const
{
	Q_RANGE_CHECK(offset, 0, nz() - 1);

	return(d[offset].value);
}

// **************************************************************************
// 
// **************************************************************************
void svec::push_nz(int idx, const mpq &v)
{
	Q_RANGE_CHECK(idx, 0, n - 1);

	d.push_back(svec_el(idx, v));
	
	if (offs_active())
		offs[idx] = nz();
}

void svec::push(int idx, const mpq &v)
{
	Q_RANGE_CHECK(idx, 0, n - 1);

	if (v.sign())
		push_nz(idx, v);
}

const mpq &svec::get(int idx) const
{
	Q_RANGE_CHECK(idx, 0, n - 1);

	int l = locate(idx);
	if (l == -1)
		return(zero);
	return(d[l].value);
}

const mpq &svec::set(int idx, const mpq &v)
{
	Q_RANGE_CHECK(idx, 0, n - 1);

	int l = locate(idx);

	if (v.sign()) {
		if (l == -1) {
			push_nz(idx, v);
			return(d[nz() - 1].value);
		}
		
		d[l].value = v;
		return(d[l].value);
	} else {
		if (l != -1)
			remove(l);
		return(zero);
	}
	
}

// **************************************************************************
// 
// **************************************************************************
const mpq &svec::operator[](int idx) const
{
	return(get(idx));
}

svec_ref svec::operator[](int idx)
{
	return(svec_ref(*this, idx));
}


// **************************************************************************
// 
// **************************************************************************
const svec &svec::operator+() const
{
	return(*this);
}

svec svec::operator-() const
{
	svec r(*this);
	
	for (int l = 0; l < nz(); l++)
		r.d[l].value.set_neg();
	
	return(r);
}


svec svec::operator+(const svec &b) const
{
	Q_MATCH(n, b.n);
	
	svec r(*this);

	r += b;
	
	return(r);
}

svec svec::operator-(const svec &b) const
{
	Q_MATCH(n, b.n);

	svec r(*this);

	r -= b;
	
	return(r);
}


// **************************************************************************
// 
// **************************************************************************
mpq svec::operator*(const svec &b) const
{
	Q_MATCH(n, b.n);

	mpq v;
	
	if (offs_active()) {
		for (int l = 0; l < b.nz(); l++) {
			int i = b.d[l].index;
			if ((i < n) && (offs[i]))
				v += d[offs[i] - 1].value * b.d[l].value;
		}
	} else if (b.offs_active()) {
		for (int l = 0; l < nz(); l++) {
			int i = d[l].index;
			if ((i < b.n) && (b.offs[i]))
				v += d[l].value * b.d[b.offs[i] - 1].value;
		}
	} else {
		dvec db;
		
		b.spread(db);
		
		v  = (*this) * db;
	}

	return(v);
}

mpq svec::operator*(const dvec &b) const
{
	Q_MATCH(n, b.size());

	mpq v;
	
	for (int l = 0; l < nz(); l++) {
		int j = d[l].index;
		
		if (j < b.size())
			v += d[l].value * b[j];
	}
	
	return(v);
}


// **************************************************************************
// 
// **************************************************************************
svec svec::operator*(const mpq &v)
{
	svec r(n);
	
	if (!v.sign())
		return(r);
	
	for (int l = 0; l < nz(); l++)
		r.push(d[l].index, d[l].value * v);
	
	return(r);
}

svec svec::operator/(const mpq &v)
{
	svec r(n);
	
	for (int l = 0; l < nz(); l++)
		r.push(d[l].index, d[l].value / v);
	
	return(r);
}

// **************************************************************************
// 
// **************************************************************************
svec &svec::operator=(const svec &b)
{
	n = b.n;
	d = b.d;
	
	if ((offs_active()) && (b.offs_active()))
		offs = b.offs;
	else
		offs_clear();
	
	return(*this);
}

svec &svec::operator+=(const svec &b)
{
	Q_MATCH(n, b.n);

	if (offs_active()) {
		for (int l = 0; l < b.nz(); l++) {
			int i = b.d[l].index;
			
			if ((i < n) && (offs[i]))
				d[offs[i] - 1].value += b.d[l].value;
			else
				push(i, b.d[l].value);
		}
	} else {
		dvec da, db;
		
		spread(da);
		b.spread(db);
		
		da += db;
		
		gather(da);
	}

	return(*this);
}

svec &svec::operator-=(const svec &b)
{
	Q_MATCH(n, b.n);

	if (offs_active()) {
		for (int l = 0; l < b.nz(); l++) {
			int i = b.d[l].index;
			
			if ((i < n) && (offs[i]))
				d[offs[i] - 1].value -= b.d[l].value;
			else
				push(i, -b.d[l].value);
		}
	} else {
		dvec da, db;
		
		spread(da);
		b.spread(db);
		
		da += db;
		
		gather(da);
	}

	return(*this);
}

svec &svec::operator*=(const mpq &v)
{
	if (!v.sign()) {
		set_zero();
		return(*this);
	}
	
	for (int l = 0; l < nz(); l++)
		d[l].value *= v;
	
	return(*this);
}

svec &svec::operator/=(const mpq &v)
{
	for (int l = 0; l < nz(); l++)
		d[l].value /= v;

	return(*this);
}

// **************************************************************************
// 
// **************************************************************************
static bool comp(const svec_el &a, const svec_el &b)
{
	return(a.index < b.index);
}

void svec::sort()
{
	std::sort(d.begin(), d.end(), comp);
}

// **************************************************************************
// 
// **************************************************************************
void svec::fdump(FILE *f, const std::string &name) const
{
	if (name.length())
		gmp_fprintf(f, "%s = ", name.c_str());
	
	for (int l = 0; l < nz(); l++) {
		gmp_fprintf(f, " %d:%Zd/%Zd",
			d[l].index,
			mpq_numref(d[l].value.v),
			mpq_denref(d[l].value.v)
			);
	}
	
	gmp_fprintf(f, "\n");
}

void svec::dump(const std::string &name) const
{
	fdump(stdout, name);
}

// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const svec &v)
{
	for (int l = 0; l < v.nz(); l++)
		os << " " << v.index(l) << ":" << v.value(l);
	
	return(os);
}


// **************************************************************************
// 
// **************************************************************************
}
