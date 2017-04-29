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
#ifndef QXX_SMAT_HPP
#define QXX_SMAT_HPP
#include <cstdio>
#include "rational.hpp"
#include "dvec.hpp"
#include "dmat.hpp"
#include "svec.hpp"


// **************************************************************************
// 
// **************************************************************************
namespace q {

// **************************************************************************
// 
// **************************************************************************
class smat {

public:
	smat();
	smat(const smat &b);
	smat(int m, int n);

	int rows() const;
	int cols() const;
	void resize(int m, int n);
	void clear();

	const mpq &get(int i, int j) const;
	void set(int i, int j, const mpq &v);
	void set_zero();
	void set_identity();

	void gather(const dmat &src);
	void spread(dmat &r) const;
	dmat dense() const;

	svec get_col(int j) const;
	void set_col(int j, const svec &v);
	const svec &get_row(int i) const;
	void set_row(int i, const svec &v);

	svec &operator[] (int i);
	const svec &operator[] (int i) const;
	svec_ref operator() (int i, int j);
	const mpq &operator() (int i, int j) const;
	
	smat &operator=(const smat &b);
	
	smat operator*(const smat &b) const;
	dvec operator*(const dvec &v) const;
	svec operator*(const svec &v) const;
	smat inv() const;
	void set_transpose(const smat &a);
	smat t() const;
	mpq det() const;

	void fdump(FILE *f, const std::string &name = "") const;
	void dump(const std::string &name = "") const;
	
private:
	int m, n;
	std::vector<svec> d;

};


// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const smat &a);


// **************************************************************************
// 
// **************************************************************************
}

#endif

