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
#ifndef QXX_TYPES_H
#define QXX_TYPES_H
#include <cstdio>
#include <iostream>
#include <gmp.h>

// **************************************************************************
// 
// **************************************************************************
namespace q {

// **************************************************************************
// 
// **************************************************************************
class mpq {

public:
	mpq();
	mpq(const mpq &b);
	mpq(long num, long den);
	mpq(int num, int den);
	mpq(long i);
	mpq(int i);
	mpq(double d);
	mpq(const char *s);
	~mpq();

	const mpq &operator+() const;
	mpq operator-() const;
	mpq operator+(const mpq &b) const;
	mpq operator-(const mpq &b) const;
	mpq operator*(const mpq &b) const;
	mpq operator/(const mpq &b) const;

	mpq &operator=(const mpq &b);
	mpq &operator=(int i);
	mpq &operator+=(const mpq &b);
	mpq &operator-=(const mpq &b);
	mpq &operator*=(const mpq &b);
	mpq &operator/=(const mpq &b);

	bool operator<(const mpq &b) const;
	bool operator>(const mpq &b) const;
	bool operator<=(const mpq &b) const;
	bool operator>=(const mpq &b) const;
	bool operator==(const mpq &b) const;
	bool operator!=(const mpq &b) const;

	void set_neg();
	void set_neg(const mpq &b);
	void set_inv();
	void set_inv(const mpq &b);
	void set_abs();
	void set_abs(const mpq &b);
	
	mpq neg() const;
	mpq inv() const;
	mpq abs() const;
	int sign() const;

	mpq num() const;
	mpq den() const;
	int bits() const;

	mpq frac() const;
	mpq floor() const;
	mpq trunc() const;
	mpq ceil() const;

	long get_long_num() const;
	long get_long_den() const;
	double get_double() const;
	mpq reduce(mpq max_den) const;

	void fdump(FILE *f, const std::string &name = "") const;
	void dump(const std::string &name = "") const;

public:
	mpq_t v;

};


// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const mpq &v);



// **************************************************************************
// 
// **************************************************************************
}

#endif


