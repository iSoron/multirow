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
#include <cassert>
#include "qxx/rational.hpp"
#include "qxx/strprintf.hpp"

// **************************************************************************
// 
// **************************************************************************
namespace q {

// **************************************************************************
// 
// **************************************************************************
mpq::mpq()
{
	mpq_init(v);
}

mpq::mpq(const mpq &b)
{
	mpq_init(v);
	mpq_set(v, b.v);
}

mpq::mpq(long num, long den)
{
	mpq_init(v);
	mpq_set_si(v, num, den);
	mpq_canonicalize(v);
}

mpq::mpq(int num, int den)
{
	mpq_init(v);
	mpq_set_si(v, num, den);
	mpq_canonicalize(v);
}

mpq::mpq(long i)
{
	mpq_init(v);
	mpq_set_si(v, i, 1);
}

mpq::mpq(int i)
{
	mpq_init(v);
	mpq_set_si(v, i, 1);
}

mpq::mpq(double d)
{
	mpq_init(v);
	mpq_set_d(v, d);
}

mpq::mpq(const char *s)
{
	mpq_init(v);
	mpq_set_str(v, s, 0);
	mpq_canonicalize(v);
}

mpq::~mpq()
{
	mpq_clear(v);
}

// **************************************************************************
// 
// **************************************************************************
const mpq &mpq::operator+() const
{
	return(*this);
}

mpq mpq::operator-() const
{
	mpq r;
	mpq_neg(r.v, v);
	return(r);
}

mpq mpq::operator+(const mpq &b) const
{
	mpq r;
	mpq_add(r.v, v, b.v);
	return(r);
}

mpq mpq::operator-(const mpq &b) const
{
	mpq r;
	mpq_sub(r.v, v, b.v);
	return(r);
}

mpq mpq::operator*(const mpq &b) const
{
	mpq r;
	mpq_mul(r.v, v, b.v);
	return(r);
}

mpq mpq::operator/(const mpq &b) const
{
	mpq r;
	mpq_div(r.v, v, b.v);
	return(r);
}


mpq &mpq::operator=(const mpq &b)
{
	mpq_set(v, b.v);
	return(*this);
}

mpq &mpq::operator=(int i)
{
	mpq_set_si(v, i, 1);
	return(*this);
}

mpq &mpq::operator+=(const mpq &b)
{
	mpq_add(v, v, b.v);
	return(*this);
}

mpq &mpq::operator-=(const mpq &b)
{
	mpq_sub(v, v, b.v);
	return(*this);
}

mpq &mpq::operator*=(const mpq &b)
{
	mpq_mul(v, v, b.v);
	return(*this);
}

mpq &mpq::operator/=(const mpq &b)
{
	mpq_div(v, v, b.v);
	return(*this);
}


bool mpq::operator<(const mpq &b) const
{
	return(mpq_cmp(v, b.v) < 0);
}

bool mpq::operator>(const mpq &b) const
{
	return(mpq_cmp(v, b.v) > 0);
}

bool mpq::operator<=(const mpq &b) const
{
	return(mpq_cmp(v, b.v) <= 0);
}

bool mpq::operator>=(const mpq &b) const
{
	return(mpq_cmp(v, b.v) >= 0);
}

bool mpq::operator==(const mpq &b) const
{
	return(mpq_equal(v, b.v));
}

bool mpq::operator!=(const mpq &b) const
{
	return(!mpq_equal(v, b.v));
}

// **************************************************************************
// 
// **************************************************************************
void mpq::set_neg()
{
	mpq_neg(v, v);
}

void mpq::set_neg(const mpq &b)
{
	mpq_neg(v, b.v);
}

void mpq::set_inv()
{
	mpq_inv(v, v);
}

void mpq::set_inv(const mpq &b)
{
	mpq_inv(v, b.v);
}

void mpq::set_abs()
{
	mpq_abs(v, v);
}

void mpq::set_abs(const mpq &b)
{
	mpq_abs(v, b.v);
}

mpq mpq::neg() const
{
	mpq r;
	mpq_neg(r.v, v);
	return(r);
}

mpq mpq::inv() const
{
	mpq r;
	mpq_inv(r.v, v);
	return(r);
}

mpq mpq::abs() const
{
	mpq r;
	mpq_abs(r.v, v);
	return(r);
}

int mpq::sign() const
{
	return(mpq_sgn(v));
}

// **************************************************************************
// 
// **************************************************************************
mpq mpq::num() const
{
	mpq r;

	mpz_set(mpq_numref(r.v), mpq_numref(v));
	mpz_set_ui(mpq_denref(r.v), 1);

	return(r);
}

mpq mpq::den() const
{
	mpq r;

	mpz_set(mpq_numref(r.v), mpq_denref(v));
	mpz_set_ui(mpq_denref(r.v), 1);

	return(r);
}

int mpq::bits() const
{
	return(
		mpz_sizeinbase(mpq_numref(v), 2)
		+ mpz_sizeinbase(mpq_denref(v), 2)
		);
}

mpq mpq::frac() const
{
	mpq r;

	mpz_fdiv_r(mpq_numref(r.v), mpq_numref(v), mpq_denref(v));
	mpz_set(mpq_denref(r.v), mpq_denref(v));

	return(r);
}

mpq mpq::floor() const
{
	mpq r;

	mpz_fdiv_q(mpq_numref(r.v), mpq_numref(v), mpq_denref(v));
	mpz_set_ui(mpq_denref(r.v), 1);

	return(r);
}

mpq mpq::trunc() const
{
	mpq r;

	mpz_tdiv_q(mpq_numref(r.v), mpq_numref(v), mpq_denref(v));
	mpz_set_ui(mpq_denref(r.v), 1);

	return(r);
}

mpq mpq::ceil() const
{
	mpq r;

	mpz_cdiv_q(mpq_numref(r.v), mpq_numref(v), mpq_denref(v));
	mpz_set_ui(mpq_denref(r.v), 1);

	return(r);
}


// **************************************************************************
// 
// **************************************************************************
long mpq::get_long_num() const
{
	assert(mpz_fits_slong_p(mpq_numref(v)));
	
	return(mpz_get_si(mpq_numref(v)));
}

long mpq::get_long_den() const
{
	assert(mpz_fits_slong_p(mpq_denref(v)));
	
	return(mpz_get_si(mpq_denref(v)));
}

double mpq::get_double() const
{
	return(mpq_get_d(v));
}

// **************************************************************************
// 
// **************************************************************************
mpq mpq::reduce(mpq max_coeff) const
{
	mpq p[3], q[3];
	mpq f, ff, r;
	
	p[0] = 0;
	q[0] = 1;
	p[1] = 1;
	q[1] = 0;
	f = *this;
	
	while (1) {
		ff = f.floor();
		r = f - ff;
		
		p[2] = ff * p[1] + p[0];
		q[2] = ff * q[1] + q[0];
		
		if ((p[2].abs() > max_coeff) || (q[2] > max_coeff))
			break;
		
		if (r.sign() == 0) {
			p[1] = p[2];
			q[1] = q[2];
			break;
		}
		
		f = r.inv();
		
		p[0] = p[1];
		q[0] = q[1];
		p[1] = p[2];
		q[1] = q[2];
	}
	

	if (q[1].sign() == 0)
		q[1] = 1;
	
	return(p[1] / q[1]);
}

// **************************************************************************
// 
// **************************************************************************
void mpq::fdump(FILE *f, const std::string &name) const
{
	mpq_t t;
	
	mpq_init(t);
	mpq_set(t, v);
	mpq_canonicalize(t);
	
	if (mpq_cmp(t, v) != 0) {
		gmp_fprintf(f, "\nBUG: not canon.: %Qd\n",
			mpq_numref(v),
			mpq_denref(v));
	}
	
	mpq_clear(t);

	if (name.length())
		gmp_fprintf(f, "%s = %Qd;\n", name.c_str(), v);
	else
		gmp_fprintf(f, " %Qd", v);
}

void mpq::dump(const std::string &name) const
{
	fdump(stdout, name);
}

// **************************************************************************
// 
// **************************************************************************
std::ostream &operator<<(std::ostream &os, const mpq &v)
{
	os << strprintf("%Qd", v.v);
	return(os);
}

// **************************************************************************
// 
// **************************************************************************
}


