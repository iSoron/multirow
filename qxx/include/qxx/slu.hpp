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
#ifndef QXX_SLU_H
#define QXX_SLU_H
#include <vector>
#include "qxx/smat.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {


// **************************************************************************
// 
// **************************************************************************
class perm {

public:
	perm();
	~perm();
	
	int size() const;
	void resize(int n);
	void clear();
	
	void id();
	void id(int n);
	
	void pivot(int k, int i);
	
public:
	std::vector<int> fwd;
	std::vector<int> bwd;

};

// **************************************************************************
// 
// **************************************************************************
class slu {

public:
	slu(const smat &a);
	~slu();
	
	dvec solve_Ax(const dvec &b) const;
	dvec solve_xA(const dvec &b) const;
	mpq det() const;

public:
	void pivot(int k);
	void factorize(const smat &a);
	dvec solve_gen(const smat &l, const smat &u,
		const perm &row, const perm &col, const dvec &b) const;
	
	int n;
	
	std::vector<int> cnz;
	perm row, col;
	smat w, lt, u;
	smat l, ut;
};

// **************************************************************************
// 
// **************************************************************************
};

#endif


