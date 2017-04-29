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
#ifndef QXX_SPARSE_VECTOR_HPP
#define QXX_SPARSE_VECTOR_HPP
#include <cstdio>
#include "rational.hpp"
#include "dvec.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {

class sparse_el;
class svec_ref;
class svec;

// **************************************************************************
// 
// **************************************************************************
class svec_ref {

public:
	svec_ref(svec &vec, int i);
	svec_ref(const svec_ref &src);
	~svec_ref();
	
	const mpq &get() const;
	
	operator mpq() const;
	
	const mpq &operator=(const svec_ref &src) const;
	const mpq &operator=(const mpq &b) const;
	const mpq &operator+=(const mpq &b) const;
	const mpq &operator-=(const mpq &b) const;
	const mpq &operator*=(const mpq &b) const;
	const mpq &operator/=(const mpq &b) const;

	const mpq &operator+() const;
	mpq operator-() const;
	mpq operator+(const mpq &b) const;
	mpq operator-(const mpq &b) const;
	mpq operator*(const mpq &b) const;
	mpq operator/(const mpq &b) const;

	bool operator<(const mpq &b) const;
	bool operator>(const mpq &b) const;
	bool operator<=(const mpq &b) const;
	bool operator>=(const mpq &b) const;
	bool operator==(const mpq &b) const;
	bool operator!=(const mpq &b) const;

public:
	svec &vector;
	int index;
};


// **************************************************************************
// 
// **************************************************************************
class svec_el {

public:
	svec_el();
	svec_el(int idx, const mpq &v);
	svec_el(const svec_el &e);
	~svec_el();

public:
	int index;
	mpq value;
};
	

// **************************************************************************
// 
// **************************************************************************
class svec {

public:
	svec();
	svec(const svec &src);
	svec(int n);
	svec(const dvec &src);
	~svec();

	void clear();
	
	int size() const;
	void resize(int n);
	int nz() const;
	
	bool offs_active() const;
	void offs_build();
	void offs_clear();

	void set_zero();

	void gather(const dvec &src);
	void spread(dvec &r) const;
	dvec dense() const;
	
	int locate(int idx) const;
	void remove(int offset);
	
	int &index(int offset);
	mpq &value(int offset);
	int index(int offset) const;
	const mpq &value(int offset) const;

	void push_nz(int idx, const mpq &v);
	void push(int idx, const mpq &v);
	
	const mpq &get(int idx) const;
	const mpq &set(int idx, const mpq &v);
	
	const mpq &operator[](int idx) const;
	svec_ref operator[](int idx);
	
	const svec &operator+() const;
	svec operator-() const;

	svec operator+(const svec &b) const;
	svec operator-(const svec &b) const;
	mpq operator*(const svec &b) const;
	mpq operator*(const dvec &b) const;
	
	svec operator*(const mpq &v);
	svec operator/(const mpq &v);
	
	svec &operator=(const svec &b);
	svec &operator+=(const svec &b);
	svec &operator-=(const svec &b);
	svec &operator*=(const mpq &v);
	svec &operator/=(const mpq &v);
	
	void sort();
	
	void fdump(FILE *f, const std::string &name = "") const;
	void dump(const std::string &name = "") const;
	
private:
	int n;
	std::vector<svec_el> d;
	std::vector<int> offs;

	static const mpq zero;
};

// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const svec &v);


// **************************************************************************
// 
// **************************************************************************
}

#endif


