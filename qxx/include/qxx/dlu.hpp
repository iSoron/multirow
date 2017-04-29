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
#ifndef QXX_DENSE_LU_HPP
#define QXX_DENSE_LU_HPP
#include "rational.hpp"
#include "dmat.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {


// **************************************************************************
// 
// **************************************************************************
class dlu {

public:
	dlu(const dmat &a);
	~dlu();

	dmat L() const;
	dmat U() const;
	dvec solve_Ax(const dvec &b) const;
	dvec solve_xA(const dvec &b) const;
	mpq det() const;

private:
	void pivot(int k);
	void factorize(const dmat &a);

	int n;
	std::vector<int> row_bwd;
	dmat g;

};

// **************************************************************************
// 
// **************************************************************************
};

#endif


