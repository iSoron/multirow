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
#ifndef QXX_DMAT_HPP
#define QXX_DMAT_HPP
#include <cstdio>
#include "rational.hpp"
#include "dvec.hpp"


// **************************************************************************
// 
// **************************************************************************
namespace q {

// **************************************************************************
// 
// **************************************************************************
class dmat {

public:
	dmat();
	dmat(const dmat &b);
	dmat(int m, int n);

	int rows() const;
	int cols() const;
	void resize(int m, int n);
	void clear();

	const mpq &get(int i, int j) const;
	void set(int i, int j, const mpq &v);
	void set(const mpq &v);
	void set_identity();

	dvec get_col(int j) const;
	void set_col(int j, const dvec &v);
	const dvec &get_row(int i) const;
	void set_row(int i, const dvec &v);

	dvec &operator[] (int i);
	const dvec &operator[] (int i) const;
	mpq &operator() (int i, int j);
	const mpq &operator() (int i, int j) const;
	
	const dmat &operator+() const;
	dmat operator-() const;
	dmat operator+(const dmat &b) const;
	dmat operator-(const dmat &b) const;
	dmat operator*(const dmat &b) const;
	dvec operator*(const dvec &v) const;
	dmat operator*(const mpq &v) const;
	dmat operator/(const mpq &v) const;
	dmat inv() const;
	dmat t() const;
	mpq det() const;

	bool operator==(const dmat &b) const;
	bool operator!=(const dmat &b) const;

	dmat &operator=(const dmat &b);
	dmat &operator=(const dvec &v);
	
	dmat &operator+=(const dmat &b);
	dmat &operator-=(const dmat &b);
	dmat &operator*=(const dmat &b);
	dmat &operator*=(const mpq &v);
	dmat &operator/=(const mpq &v);
	
	void fdump(FILE *f, const std::string &name = "") const;
	void dump(const std::string &name = "") const;
	
private:
	int m, n;
	std::vector<dvec> d;

};

// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const dmat &a);


// **************************************************************************
// 
// **************************************************************************
}

#endif

