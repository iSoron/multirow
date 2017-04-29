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
#ifndef QXX_DENSE_VECTOR_HPP
#define QXX_DENSE_VECTOR_HPP
#include <cstdio>
#include <vector>
#include "rational.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {

class svec;

// **************************************************************************
// 
// **************************************************************************
class dvec {

public:
	dvec();
	dvec(const dvec &b);
	dvec(int n);
	dvec(const svec &b);
	
	int size() const;
	void resize(int n);
	void clear();
	mpq &operator[] (int i);
	const mpq &operator[] (int i) const;
	
	void set(const mpq &v);
	void set(int i0, int i1, const mpq &v);
	void assign(int n, const mpq &v);
	
	dvec operator+(const dvec &b) const;
	dvec operator-(const dvec &b) const;
	mpq operator*(const dvec &b) const;
	mpq operator*(const svec &b) const;
	
	const dvec &operator+() const;
	dvec operator-() const;
	dvec operator+(const mpq &v) const;
	dvec operator-(const mpq &v) const;
	dvec operator*(const mpq &v) const;
	dvec operator/(const mpq &v) const;
	
	bool operator==(const dvec &b) const;
	bool operator!=(const dvec &b) const;
	
	dvec &operator=(const dvec &b);
	dvec &operator=(const svec &b);

	dvec &operator+=(const dvec &b);
	dvec &operator-=(const dvec &b);
	dvec &operator+=(const mpq &v);
	dvec &operator-=(const mpq &v);
	dvec &operator*=(const mpq &v);
	dvec &operator/=(const mpq &v);
	
	void fdump(FILE *f, const std::string &name = "") const;
	void dump(const std::string &name = "") const;

private:
	std::vector<mpq> d;

};


// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const dvec &v);

// **************************************************************************
// 
// **************************************************************************
}

#endif


